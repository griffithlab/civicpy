from civicpy import civic
from csv import DictWriter
import datetime
from civicpy.__version__ import __version__


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

    SIMPLE_VARIANT_TYPES = [
        'SO:0001583',
        'SO:0001910',
        'SO:0001587',
        'SO:0001822',
        'SO:0001818',
        'SO:0001589',
    ]

    VALID_VARIANTS = dict()

    def __init__(self, f, version=4.2):
        self._f = f
        assert version in [4.2]  # Supported VCF versions
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
            # Requires the variant type to be specified by CIViC
            types = variant.types
            if not types:
                self.VALID_VARIANTS[variant] = False
                return
            # Requires the variant type(s) to be SNV/indels
            if not all([x.so_id in self.SIMPLE_VARIANT_TYPES for x in types]):
                self.VALID_VARIANTS[variant] = False
                return
            # Currently requires exactly one coordinate set

            # Passes all above validations
            self.VALID_VARIANTS[variant] = True
        elif not valid:
            return
