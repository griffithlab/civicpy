from civicpy import civic
from csv import DictWriter
import datetime
from civicpy.__version__ import __version__
import requests


class VCFWriter(DictWriter):
    """
    :param filehandle f: A filehandle for the VCF output file
    """
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

    def __init__(self, f, version=4.2):
        self._f = f
        assert version in VCFWriter.SUPPORTED_VERSIONS  # Supported VCF versions
        self.version = version
        super().__init__(f, delimiter='\t', fieldnames=self.HEADER, restval='.', lineterminator='\n')
        self.meta_info_fields = []
        self.variant_records = set()

    def writeheader(self):
        """
        Writes the header lines to the VCF file.
        """
        # write meta lines
        self._write_meta_file_lines()
        self._write_meta_info_lines()
        # write header line
        super().writeheader()

    def addrecord(self, civic_record):
        """
        Takes either a :class:`civic.Evidence`, :class:`civic.Assertion`, :class:`civic.Variant`, or :class:`civic.Gene` object
        and adds all :class:`civic.Variant` objects associated with it to the VCFWriter object for processing and writing to the VCF. 

        :param civic.CivicRecord civic_record: Either a :class:`civic.Evidence`, :class:`civic.Assertion`, :class:`civic.Variant`, or :class:`civic.Gene` object
        """
        if isinstance(civic_record, civic.Evidence) or isinstance(civic_record, civic.Assertion):
            if civic_record.variant.is_valid_for_vcf(emit_warnings=True):
                self._add_variant_record(civic_record.variant)
        elif isinstance(civic_record, civic.Gene):
            for variant in civic_record.variants:
                if civic_record.variant.is_valid_for_vcf(emit_warnings=True):
                    self._add_variant_record(variant)
        elif isinstance(civic_record, civic.Variant):
            if civic_record.is_valid_for_vcf(emit_warnings=True):
                self._add_variant_record(civic_record)
        else:
            raise ValueError('Expected a CIViC Variant, Assertion or Evidence record.')

    def addrecords(self, civic_records):
        """
        Takes multiple :class:`civic.Evidence`, :class:`civic.Assertion`, :class:`civic.Variant`, and/or :class:`civic.Gene` objects
        and adds all :class:`civic.Variant` objects associated with them to the VCFWriter object for processing and writing to the VCF.
        ``civic_records`` can contain a mix of these object types.

        :param list civic_records: A list of a :class:`civic.Evidence`, :class:`civic.Assertion`, :class:`civic.Variant`, and/or :class:`civic.Gene` objects
        """
        for record in civic_records:
            self.addrecord(record)

    def writerecords(self, with_header=True):
        """
        Takes all variant objects saved to the VCFWriter object, processes them, and outputs them to the VCF file

        :param bool with_header: Indicates weather or not the VCF header lines should be written as part of this function call.
        """
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

        ensembl_server = "https://grch37.rest.ensembl.org"

        # write them
        rows = []
        for variant in sorted_records:
            if variant.coordinates.reference_build != 'GRCh37':
                continue
            if variant.is_insertion:
                if not variant.coordinates.representative_transcript:
                    continue
                else:
                    start = variant.coordinates.start
                    ext = "/sequence/region/human/{}:{}-{}".format(variant.coordinates.chromosome, start, start)
                    r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
                    if variant.coordinates.reference_bases == None or variant.coordinates.reference_bases == '-' or variant.coordinates.reference_bases == '':
                        ref = r.text
                    else:
                        ref = "{}{}".format(r.text, variant.coordinates.reference_bases)
                    alt = "{}{}".format(r.text, variant.coordinates.variant_bases)
                    csq_alt = variant.coordinates.variant_bases
            elif variant.is_deletion:
                if not variant.coordinates.representative_transcript:
                    continue
                else:
                    start = variant.coordinates.start - 1
                    ext = "/sequence/region/human/{}:{}-{}".format(variant.coordinates.chromosome, start, start)
                    r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
                    ref = "{}{}".format(r.text, variant.coordinates.reference_bases)
                    if variant.coordinates.variant_bases == None or variant.coordinates.variant_bases == '-' or variant.coordinates.variant_bases == '':
                        alt = r.text
                    else:
                        alt = "{}{}".format(r.text, variant.coordinates.variant_bases)
                    csq_alt = "-"
            else:
                start = variant.coordinates.start
                ref = variant.coordinates.reference_bases
                alt = variant.coordinates.variant_bases
                csq_alt = alt
            out_dict = {
                '#CHROM': variant.coordinates.chromosome,
                'POS':    str(start),
                'ID':     variant.id,
                'REF':    ref,
                'ALT':    alt,
            }

            info_dict = {
                'GN': variant.gene.name,
                'VT': variant.name,
            }
            csq = []
            if variant.coordinates.representative_transcript:
                hgvs_cs = [e for e in variant.hgvs_expressions if (':c.' in e) and (variant.coordinates.representative_transcript in e)]
                hgvs_ps = [e for e in variant.hgvs_expressions if (':p.' in e) and (variant.coordinates.representative_transcript in e)]
            hgvs_c = hgvs_cs[0] if len(hgvs_cs) == 1 else ''
            hgvs_p = hgvs_ps[0] if len(hgvs_ps) == 1 else ''
            for evidence in variant.evidence:
                special_character_table = str.maketrans(VCFWriter.SPECIAL_CHARACTERS)
                csq.append('|'.join([
                    csq_alt,
                    '&'.join(map(lambda t: t.name, variant.variant_types)),
                    variant.gene.name,
                    str(variant.gene.entrez_id),
                    'transcript',
                    str(variant.coordinates.representative_transcript),
                    hgvs_c,
                    hgvs_p,
                    variant.name,
                    str(variant.id),
                    '&'.join(map(lambda a: a.translate(special_character_table), variant.variant_aliases)),
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
                    csq_alt,
                    '&'.join(map(lambda t: t.name, variant.variant_types)),
                    variant.gene.name,
                    str(variant.gene.entrez_id),
                    'transcript',
                    str(variant.coordinates.representative_transcript),
                    hgvs_c,
                    hgvs_p,
                    variant.name,
                    str(variant.id),
                    '&'.join(map(lambda a: a.translate(special_character_table), variant.variant_aliases)),
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
            rows.append(out_dict)
        return rows

    def _write_meta_file_lines(self):
        self._f.write(f'##fileformat=VCFv{self.version}\n')
        self._f.write('##fileDate={}\n'.format(
            datetime.date.today().strftime('%Y%m%d')
        ))
        self._f.write('##reference=ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz\n')
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
            'Consequence',
            'SYMBOL',
            'Entrez Gene ID',
            'Feature_type',
            'Feature',
            'HGVSc',
            'HGVSp',
            'CIViC Variant Name',
            'CIViC Variant ID',
            'CIViC Variant Aliases',
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
        s = [f'ID={id_},Number={number},Type={type_},Description="{description}"']
        s.extend([f'{k}={v}' for k, v in kwargs])
        out = ','.join(s)
        self._f.write(f'##INFO=<{out}>\n')

    def _add_variant_record(self, variant_record):
        self.variant_records.add(variant_record)
