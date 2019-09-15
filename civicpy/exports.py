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
    SPECIAL_CHARACTERS = {
        " ": r'\x20',
        "!": r'\x21',
        '"': r'\x22',
        '#': r'\x23',
        '$': r'\x24',
        '%': r'\x25',
        '&': r'\x26',
        "'": r'\x27',
        '(': r'\x28',
        ')': r'\x29',
        '*': r'\x2A',
        '+': r'\x2B',
        ',': r'\x2C',
        '-': r'\x2D',
        '.': r'\x2E',
        '/': r'\x2F',
        ':': r'\x3A',
        ';': r'\x3B',
        '<': r'\x3C',
        '=': r'\x3D',
        '>': r'\x3E',
        '?': r'\x3F',
        '@': r'\x40',
        '[': r'\x5B',
        '\\': r'\x5C',
        ']': r'\x5D',
        '^': r'\x5E',
        '_': r'\x5F',
        '`': r'\x60',
        '{': r'\x7B',
        '|': r'\x7C',
        '}': r'\x7D',
        '~': r'\x7E',
    }

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
        self.variant_records = set()

    def writeheader(self):
        # write meta lines
        self._write_meta_file_lines()
        self._write_meta_info_lines()
        # write header line
        super().writeheader()

    def addrecord(self, civic_record):
        if isinstance(civic_record, civic.Evidence) or isinstance(civic_record, civic.Assertion):
            self._validate_variant(civic_record.variant)
            self._add_variant_record(civic_record.variant)
        elif isinstance(civic_record, civic.Gene):
            for variant in civic_record.variants:
                self._validate_variant(variant)
                self._add_variant_record(variant)
        elif isinstance(civic_record, civic.Variant):
            self._validate_variant(civic_record)
            self._add_variant_record(civic_record)
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
        sorted_records = list(self.variant_records)
        sorted_records.sort(key=lambda x: int(x.coordinates.stop))
        sorted_records.sort(key=lambda x: int(x.coordinates.start))
        int_chromosomes = [i for i in sorted_records if i.coordinates.chromosome.isdigit()]
        string_chromosomes = [i for i in sorted_records if not i.coordinates.chromosome.isdigit()]
        int_chromosomes.sort(key=lambda x: int(x.coordinates.chromosome))
        string_chromosomes.sort(key=lambda x: x.coordinates.chromosome)
        sorted_records = int_chromosomes + string_chromosomes

        # write them
        for variant in sorted_records:
            out_dict = {
                '#CHROM': variant.coordinates.chromosome,
                'POS':    variant.coordinates.start,
                'ID':     variant.id,
                'REF':    variant.coordinates.reference_bases,
                'ALT':    variant.coordinates.variant_bases
            }
            assert all([c.upper() in ['A', 'C', 'G', 'T', 'N'] for c in out_dict['REF']])
            assert all([c.upper() in ['A', 'C', 'G', 'T', 'N', '*', '-'] for c in out_dict['ALT']]), \
                f'observed incompatible alt allele in {evidence.variant}'

            info_dict = {
                'GN': variant.gene.name,
                'VT': variant.name,
            }
            csq = []
            for evidence in variant.evidence:
                special_character_table = str.maketrans(VCFWriter.SPECIAL_CHARACTERS)
                csq.append('|'.join([
                    out_dict['ALT'],
                    variant.gene.name,
                    str(variant.gene.entrez_id),
                    str(variant.coordinates.representative_transcript),
                    variant.name,
                    str(variant.id),
                    '&'.join(map(lambda a: a.translate(special_character_table), variant.variant_aliases)),
                    '&'.join(map(lambda t: t.name, variant.variant_types)),
                    '&'.join(map(lambda e: e.translate(special_character_table), variant.hgvs_expressions)),
                    str(variant.allele_registry_id),
                    '&'.join(variant.clinvar_entries),
                    str(variant.civic_actionability_score),
                    "evidence",
                    str(evidence.id),
                    "https://civicdb.org/links/evidence/{}".format(evidence.id),
                    "{} ({})".format(evidence.source.citation_id, evidence.source.source_type),
                    str(evidence.variant_origin),
                    evidence.status,
                ]))
            for assertion in variant.assertions:
                csq.append('|'.join([
                    out_dict['ALT'],
                    variant.gene.name,
                    str(variant.gene.entrez_id),
                    str(variant.coordinates.representative_transcript),
                    variant.name,
                    str(variant.id),
                    '&'.join(map(lambda a: a.translate(special_character_table), variant.variant_aliases)),
                    '&'.join(map(lambda t: t.name, variant.variant_types)),
                    '&'.join(map(lambda e: e.translate(special_character_table), variant.hgvs_expressions)),
                    str(variant.allele_registry_id),
                    '&'.join(variant.clinvar_entries),
                    str(variant.civic_actionability_score),
                    "assertion",
                    str(assertion.id),
                    "https://civicdb.org/links/assertion/{}".format(assertion.id),
                    "",
                    str(assertion.variant_origin),
                    assertion.status,
                ]))
            info_dict['CSQ'] = ','.join(csq)

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
        # CSQ
        self._write_meta_info_line('CSQ', '.', 'String', 'Consequence annotations from CIViC. Format: {}'.format('|'.join([
            'Allele',
            'HGNC Gene Symbol',
            'Entrez Gene ID',
            'CIViC Representative Transcript',
            'CIViC Variant Name',
            'CIViC Variant ID',
            'CIViC Variant Aliases',
            'CIViC Variant Type',
            'CIViC HGVS',
            'Allele Registry ID',
            'ClinVar IDs',
            'CIViC Variant Evidence Score',
            'CIViC Entity Type',
            'CIViC Entity ID',
            'CIViC Entity URL',
            'CIViC Entity Source',
            'CIViC Entity Variant Origin',
            'CIViC Entity Status',
        ])))

    def _write_meta_info_line(self, id_, number, type_, description, **kwargs):
        assert id_ not in self.meta_info_fields
        assert id_ not in self.VCF_RESERVED_FIELDS
        self.meta_info_fields.append(id_)
        s = [f'ID={id_},Number={number},Type={type_},Description={description}']
        s.extend([f'{k}={v}' for k, v in kwargs])
        out = ','.join(s)
        self._f.write(f'##INFO=<{out}>\n')

    def _validate_variant(self, variant):
        valid = self.VALID_VARIANTS.get(variant, None)
        if valid is None:
            valid = self._validate_sequence_variant(variant)
        if not valid:
            logging.info(f'{variant} is invalid for VCF.')
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

    def _add_variant_record(self, variant_record):
        self.variant_records.add(variant_record)
