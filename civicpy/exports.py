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

    SO_READER = SequenceOntologyReader()

    def __init__(self, f, version=4.2):
        self._f = f
        assert version in VCFWriter.SUPPORTED_VERSIONS  # Supported VCF versions
        self.version = version
        super().__init__(f, delimiter='\t', fieldnames=self.HEADER)
        self.meta_info_fields = []

    def writeheader(self):
        # write meta lines
        self._write_meta_file_lines()
        self._write_meta_info_lines()
        # write header line
        super().writeheader()

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
        # Evidence statement
        self._write_meta_info_line('ES', 1, 'String', 'CIViC Evidence Statement')
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
        self.meta_info_fields.append(id_)
        s = [f'ID={id_},Number={number},Type={type_},Description={description}']
        s.extend([f'{k}={v}' for k, v in kwargs])
        out = ','.join(s)
        self._f.write(f'##INFO=<{out}>\n')

    def writerow(self, civic_record):
        if isinstance(civic_record, civic.Variant) or isinstance(civic_record, civic.Assertion):
            for evidence in civic_record.evidence:
                self.writerow(evidence)
        elif isinstance(civic_record, civic.Evidence):
            self._validate_and_write_evidence_record(civic_record)

    def _validate_and_write_evidence_record(self, record):
        variant = record.variant
        valid = self.VALID_VARIANTS.get(variant, None)
        if valid is None:
            # valid = self._validate_structural_variant(variant)
            valid = self._validate_coordinates(variant)
        if not valid:
            logging.info(f'Skipping {record}, invalid variant {variant}.')
            return

        raise NotImplementedError  # TODO: call write evidence method here

    def _validate_structural_variant(self, variant):
        # Requires all types to have SO_IDs
        types = variant.types
        for variant_type in types:
            if not variant_type.so_id.startswith('SO:'):
                return self._cache_variant_validation(variant, False)

        # Filter types if multiple direct lineage to most specific, remove non-structural types
        type_len = len(types)
        simplified_types = list()
        for i in range(type_len):
            remove = False
            try:
                structural = self.SO_READER.same_or_has_ancestor(types[i].so_id, 'SO:0001537')
            except networkx.NetworkXError as e:
                logging.warning(f'Error for variant {variant}: {e.args[0]}')
                return self._cache_variant_validation(variant, False)
            if not structural:
                continue
            for j in range(type_len):
                if i == j:
                    continue
                if types[i] == types[j] and i > j:
                    remove = True
                elif self.SO_READER.same_or_has_descendant(types[i].so_id, types[j].so_id):
                    remove = True
            if remove:
                logging.warning(f'Simplifying redundant type {types[i].name} in {variant}')
                continue
            simplified_types.append(types[i])
        types = simplified_types

        # Requires at least one variant type (other than filtered types above) to be specified by CIViC
        return self._cache_variant_validation(variant, bool(types))

    def _validate_coordinates(self, variant):
        # Requires exactly one coordinate set with ref and alt
        coordinates = variant.coordinates
        valid = all([
            coordinates.chromosome,
            coordinates.start,
            coordinates.stop,
            coordinates.reference_bases,
            coordinates.variant_bases,
            not coordinates.chromosome2,
            not coordinates.start2,
            not coordinates.stop2
        ])
        return self._cache_variant_validation(variant, valid)

    def _cache_variant_validation(self, variant, result):
        self.VALID_VARIANTS[variant] = result
        return result