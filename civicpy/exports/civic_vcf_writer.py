from civicpy import civic
import datetime
from civicpy.__version__ import __version__
import vcfpy


class CivicVcfWriter():
    """
    :param path filepath: A file path for the VCF output file
    :param [civic.exports.CivicVcfRecord] vcf_records: A list of :class:`civic.exports.CivicVcfRecord` objects to write to the output file
    :param string version: VCF version
    """

    CSQ_DESCRIPTION = 'Consequence annotations from CIViC. Format: {}'.format('|'.join([
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
        'CIViC Variant URL',
        'CIViC Molecular Profile Name',
        'CIViC Molecular Profile ID',
        'CIViC Molecular Profile Aliases',
        'CIViC Molecular Profile URL',
        'CIViC HGVS',
        'Allele Registry ID',
        'ClinVar IDs',
        'CIViC Molecular Profile Score',
        'CIViC Entity Type',
        'CIViC Entity ID',
        'CIViC Entity URL',
        'CIViC Entity Source',
        'CIViC Entity Variant Origin',
        'CIViC Entity Status',
        'CIViC Entity Significance',
        'CIViC Entity Direction',
        'CIViC Entity Disease',
        'CIViC Entity Therapies',
        'CIViC Entity Therapy Interaction Type',
        'CIViC Evidence Phenotypes',
        'CIViC Evidence Level',
        'CIViC Evidence Rating',
        'CIViC Assertion ACMG Codes',
        'CIViC Assertion AMP Category',
        'CIViC Assertion NCCN Guideline',
        'CIViC Assertion Regulatory Approval',
        'CIViC Assertion FDA Companion Test',
    ]))

    SUPPORTED_VERSIONS = [4.2]

    def __init__(self, filepath, vcf_records, version=4.2):
        assert version in self.SUPPORTED_VERSIONS  # Supported VCF versions
        self.version = version
        self.create_header()
        vcf_writer = vcfpy.Writer.from_path(filepath, self.header)
        sorted_vcf_records = self.sort_vcf_records(vcf_records)
        for vcf_record in sorted_vcf_records:
            vcf_writer.write_record(vcf_record)
        vcf_writer.close()

    def create_header(self):
        self.header = vcfpy.Header(samples=vcfpy.SamplesInfos([]))
        self.add_metadata_header_lines()
        self.add_info_header_lines()

    def add_metadata_header_lines(self):
        self.header.add_line(vcfpy.HeaderLine('fileformat', 'VCFv{}'.format(self.version)))
        self.header.add_line(vcfpy.HeaderLine('fileDate', datetime.date.today().strftime('%Y%m%d')))
        self.header.add_line(vcfpy.HeaderLine('reference', 'ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz'))
        self.header.add_line(vcfpy.HeaderLine('source', 'CIViCpy_v{}'.format(__version__)))
        self.header.add_line(vcfpy.HeaderLine('aboutURL', 'https://civicdb.org/help/evidence/overview'))

    def add_info_header_lines(self):
        # Gene
        self.header.add_info_line(vcfpy.OrderedDict([
            ('ID', 'GN'),
            ('Number', '1'),
            ('Type', 'String'),
            ('Description', 'HGNC Gene Symbol'),
        ]))
        # Variant
        self.header.add_info_line(vcfpy.OrderedDict([
            ('ID', 'VT'),
            ('Number', '1'),
            ('Type', 'String'),
            ('Description', 'CIViC Variant Name'),
        ]))
        # CSQ
        self.header.add_info_line(vcfpy.OrderedDict([
            ('ID', 'CSQ'),
            ('Number', '.'),
            ('Type', 'String'),
            ('Description', self.CSQ_DESCRIPTION),
        ]))

    def sort_vcf_records(self, vcf_records):
        vcf_records.sort(key=lambda x: int(x.POS))
        int_chromosomes = [i for i in vcf_records if i.CHROM.isdigit()]
        string_chromosomes = [i for i in vcf_records if not i.CHROM.isdigit()]
        int_chromosomes.sort(key=lambda x: int(x.CHROM))
        string_chromosomes.sort(key=lambda x: x.CHROM)
        sorted_records = int_chromosomes + string_chromosomes
        return sorted_records
