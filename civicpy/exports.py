from civicpy import civic
from csv import DictWriter
import datetime
from civicpy.__version__ import __version__
import obonet, networkx
import logging


class SequenceOntologyReader():

    def __init__(self, url='https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so.obo'):
        self.url = url
        self.graph = obonet.read_obo(url)
        assert networkx.is_directed_acyclic_graph(self.graph)
        self.ancestor_cache = dict()

    def get_child_ids(self, id_):
        return [x for x in self.graph.predecessors(id_)]  # I know. It looks wrong, but its not.

    def get_parent_ids(self, id_):
        return [x for x in self.graph.successors(id_)]    # I know. It looks wrong, but its not.

    def same_or_has_ancestor(self, id_, ancestor_id):
        if id_ == ancestor_id:
            return True
        key = (id_, ancestor_id)
        cached = self.ancestor_cache.get(key, None)
        if cached is not None:
            return cached
        parents = self.get_parent_ids(id_)
        if not parents:
            self.ancestor_cache[key] = False
            return False
        result = any([self.same_or_has_ancestor(parent, ancestor_id) for parent in parents])
        self.ancestor_cache[key] = result
        return result

    def same_or_has_descendant(self, id_, descendant_id):
        if id_ == descendant_id:
            return True
        key = (id_, descendant_id)
        cached = self.ancestor_cache.get(key, None)
        if cached is not None:
            return cached
        children = self.get_child_ids(id_)
        if not children:
            self.ancestor_cache[key] = False
            return False
        result = any([self.same_or_has_descendant(child, descendant_id) for child in children])
        self.ancestor_cache[key] = result
        return result


class VCFWriter(DictWriter):

    HEADER = [
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    ]

    VALID_VARIANTS = dict()

    SUPPORTED_VERSIONS = [4.2]

    VCF_RESERVED_FIELDS = {
        'AA',
        'AC',
        'AF',
        'AN',
        'BQ',
        'CIGAR',
        'DB',
        'DP',
        'END',
        'H2',
        'H3',
        'MQ',
        'MQ0',
        'NS',
        'SB',
        'SOMATIC',
        'VALIDATED',
        '1000G'
    }

    SO_READER = SequenceOntologyReader()

    def __init__(self, f, version=4.2):
        self._f = f
        assert version in VCFWriter.SUPPORTED_VERSIONS  # Supported VCF versions
        self.version = version
        super().__init__(f, delimiter='\t', fieldnames=self.HEADER, restval='.', lineterminator='\n')
        self.meta_info_fields = []
        self.evidence_records = set()

    def writeheader(self):
        # write meta lines
        self._write_meta_file_lines()
        self._write_meta_info_lines()
        # write header line
        super().writeheader()

    def addrecord(self, civic_record):
        if isinstance(civic_record, civic.Variant) or isinstance(civic_record, civic.Assertion):
            for evidence in civic_record.evidence:
                self.addrecord(evidence)
        elif isinstance(civic_record, civic.Evidence):
            valid = self._validate_evidence_record(civic_record)
            if valid:
                self._add_evidence_record(civic_record)
        else:
            raise ValueError('Expected a CIViC Variant, Assertion or Evidence record.')

    def addrecords(self, civic_records):
        for record in civic_records:
            self.addrecord(record)

    def writerecords(self, with_header=True):
        # write header
        if with_header:
            self.writeheader()

        # sort records
        sorted_records = list(self.evidence_records)
        sorted_records.sort(key=lambda x: x.id)
        sorted_records.sort(key=lambda x: x.variant.coordinates.stop)
        sorted_records.sort(key=lambda x: x.variant.coordinates.start)
        sorted_records.sort(key=lambda x: x.variant.coordinates.chromosome)

        # write them
        for evidence in sorted_records:
            out_dict = {
                '#CHROM': evidence.variant.coordinates.chromosome,
                'POS':    evidence.variant.coordinates.start,
                'ID':     f'EID{evidence.id}',
                'REF':    evidence.variant.coordinates.reference_bases,
                'ALT':    evidence.variant.coordinates.variant_bases
            }
            assert all([c.upper() in ['A', 'C', 'G', 'T', 'N'] for c in out_dict['REF']])
            assert all([c.upper() in ['A', 'C', 'G', 'T', 'N', '*'] for c in out_dict['ALT']]), \
                f'observed incompatible alt allele in {evidence.variant}'

            info_dict = {
                'GN': evidence.variant.gene.name,
                'VT': evidence.variant.name,
                # 'ES': evidence.description,
                'EL': evidence.evidence_level,
                'ET': evidence.evidence_type,
                'ED': evidence.evidence_direction,
                'CS': evidence.clinical_significance,
                'VO': evidence.variant_origin,
                'DG': ','.join([d.name for d in evidence.drugs]),
                'DI': evidence.drug_interaction_type,
                'PM': evidence.source.pubmed_id,
                'TR': evidence.rating,
                'EU': evidence.site_link
            }
            doid = evidence.disease.doid
            if doid is None:
                doid = '.'
            else:
                doid = f'DOID{doid}'
            info_dict['DS'] = ','.join([evidence.disease.name, doid])

            out = list()
            for field in self.meta_info_fields:
                v = info_dict[field]
                if isinstance(v, str):
                    v = v.replace(' ', '_')
                    assert ';' not in v
                    assert '=' not in v
                if v:
                    out.append(f'{field}={v}')
            out_dict['INFO'] = ';'.join(out)

            super().writerow(out_dict)

    def _write_meta_file_lines(self):
        self._f.write(f'##fileformat=VCFv{self.version}\n')
        self._f.write('##fileDate={}\n'.format(
            datetime.date.today().strftime('%Y%m%d')
        ))
        self._f.write(f'##source=CIViCpy_v{__version__}\n')
        self._f.write(f'##aboutURL=https://civicdb.org/help/evidence/overview\n')

    def _write_meta_info_lines(self):
        # Gene
        self._write_meta_info_line('GN', 1, 'String', 'HGNC Gene Symbol')
        # Variant
        self._write_meta_info_line('VT', 1, 'String', 'CIViC Variant Name')
        # # Evidence statement
        # self._write_meta_info_line('ES', 1, 'String', 'CIViC Evidence Statement')
        # Evidence level
        self._write_meta_info_line('EL', 1, 'Character', 'CIViC Evidence Level')
        # Evidence type
        self._write_meta_info_line('ET', 1, 'String', 'CIViC Evidence Type')
        # Evidence direction
        self._write_meta_info_line('ED', 1, 'String', 'CIViC Evidence Direction')
        # Clinical significance
        self._write_meta_info_line('CS', 1, 'String', 'CIViC Clinical Significance')
        # Variant origin
        self._write_meta_info_line('VO', 1, 'String', 'CIViC Variant Origin')
        # Disease
        self._write_meta_info_line('DS', 2, 'String', 'Disease Name, DOID')
        # Drug
        self._write_meta_info_line('DG', '.', 'String', 'Drugs')
        # Drug interaction type
        self._write_meta_info_line('DI', 1, 'String', 'Interaction Between Drugs')
        # Pubmed ID
        self._write_meta_info_line('PM', 1, 'String', 'Pubmed Identifier')
        # Evidence trust rating
        self._write_meta_info_line('TR', 1, 'Integer', 'Trust rating (stars)')
        # Evidence URL
        self._write_meta_info_line('EU', 1, 'String', 'Evidence Record URL')

    def _write_meta_info_line(self, id_, number, type_, description, **kwargs):
        assert id_ not in self.meta_info_fields
        assert id_ not in self.VCF_RESERVED_FIELDS
        self.meta_info_fields.append(id_)
        s = [f'ID={id_},Number={number},Type={type_},Description={description}']
        s.extend([f'{k}={v}' for k, v in kwargs])
        out = ','.join(s)
        self._f.write(f'##INFO=<{out}>\n')

    def _validate_evidence_record(self, record):
        variant = record.variant
        valid = self.VALID_VARIANTS.get(variant, None)
        if valid is None:
            valid = self._validate_sequence_variant(variant)
        if not valid:
            logging.info(f'{record} has invalid VCF variant {variant}.')
        return valid

    def _validate_sequence_variant(self, variant):
        # Requires all types to have SO_IDs
        types = variant.types
        for variant_type in types:
            if not variant_type.so_id.startswith('SO:'):
                return self._cache_variant_validation(variant, False)

        # Filter types if multiple direct lineage to most specific, remove non-variant types
        type_len = len(types)
        simplified_types = list()
        for i in range(type_len):
            remove = False
            try:
                sequence_variant = self.SO_READER.same_or_has_ancestor(types[i].so_id, 'SO:0001060')
            except networkx.NetworkXError as e:
                logging.warning(f'Error for variant {variant}: {e.args[0]}')
                return self._cache_variant_validation(variant, False)
            if not sequence_variant:
                continue
            for j in range(type_len):
                if i == j:
                    continue
                if types[i].id == types[j].id and i > j:
                    remove = True
                elif self.SO_READER.same_or_has_descendant(types[i].so_id, types[j].so_id):
                    remove = True
            if remove:
                logging.warning(f'Simplifying redundant type {types[i].name} in {variant}')
                continue
            simplified_types.append(types[i])
        types = simplified_types

        valid = self._validate_coordinates(variant, types)

        # Requires at least one variant type (other than filtered types above) to be specified by CIViC
        return self._cache_variant_validation(variant, valid)

    def _validate_coordinates(self, variant, types):
        # If multiple types, requires exactly one to be structural type.
        if not types:
            return False

        if len(types) > 1:
            structural_types = [t for t in types if self.SO_READER.same_or_has_ancestor(t.so_id, 'SO:0001537')]
            if len(structural_types) == 1:
                types = structural_types
            elif len(structural_types) > 1:
                logging.warning(f'Variant {variant} has multiple structural types. Skipping.')
                return False
            else:
                logging.warning(f'Variant {variant} has multiple types, none structural. Skipping.')
                return False

        variant_type = types[0]

        # If type is a transcript variant, requires exactly one coordinate set with ref and alt
        if self.SO_READER.same_or_has_ancestor(variant_type.so_id, 'SO:0001576'):
            coordinates = variant.coordinates
            valid_array = [bool(x) for x in [
                coordinates.chromosome,
                coordinates.start,
                coordinates.stop,
                coordinates.reference_bases,
                coordinates.variant_bases,
                not coordinates.chromosome2,
                not coordinates.start2,
                not coordinates.stop2
            ]]
            valid = all(valid_array) \
                    and all([c.upper() in ['A', 'C', 'G', 'T'] for c in coordinates.variant_bases]) \
                    and all([c.upper() in ['A', 'C', 'G', 'T'] for c in coordinates.reference_bases])
            if not valid:
                if sum(valid_array[:5]) == 0:
                    # Nothing to do here. No inference is to be performed, and no coordinates are provided.
                    logging.warning(f'Variant {variant} has a structural type but no coordinates. Skipping.')
                elif sum(valid_array[:3]) + sum(valid_array[-3:]) == 6:
                    if sum(valid_array[3:5]) == 0:
                        # Here, neither ref nor alt is specified, as in ambiguous mutations for an amino acid.
                        logging.warning(f'Variant {variant} has a structural type but no ref or alt. Skipping.')
                    elif self.SO_READER.same_or_has_ancestor(variant_type.so_id, 'SO:0001589') or \
                        self.SO_READER.same_or_has_ancestor(variant_type.so_id, 'SO:0001820') or \
                        self.SO_READER.same_or_has_ancestor(variant_type.so_id, 'SO:0001587') or \
                        self.SO_READER.same_or_has_ancestor(variant_type.so_id, 'SO:0002012'):
                        # Here, one of ref or alt is specified, and is of a compatible variant type for an indel
                        # These are allowed.
                        valid = True
                    else:
                        raise ValueError(f'Unexpected type ({variant_type.name}) for variant ( {variant.site_link} ).')
                else:
                    raise ValueError(f'Unexpected coordinates for ( {variant.site_link} ).')
            return self._cache_variant_validation(variant, valid)
        else:
            raise NotImplementedError(f'No logic to handle {variant_type.name} {variant}')
            # TODO: handle non-transcript variants here. Currently aren't any that meet other criteria.

    def _cache_variant_validation(self, variant, result):
        self.VALID_VARIANTS[variant] = result
        return result

    def _add_evidence_record(self, evidence_record):
        self.evidence_records.add(evidence_record)
