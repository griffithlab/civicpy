import pytest
from civicpy import civic
from civicpy.exports import VCFWriter
import io


@pytest.fixture(scope='module')
def vcf_stream():
    return io.StringIO()


@pytest.fixture(scope='function')
def vcf_writer(vcf_stream):
    return VCFWriter(vcf_stream)

# snv


@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)

# simple insertion


@pytest.fixture(scope="module")
def a56fs():
    return civic.get_variant_by_id(1785)

# simple deletion


@pytest.fixture(scope="module")
def v273fs():
    return civic.get_variant_by_id(762)

# complex insertion


@pytest.fixture(scope="module")
def v2444fs():
    return civic.get_variant_by_id(137)

# complex deletion


@pytest.fixture(scope="module")
def l158fs():
    return civic.get_variant_by_id(2137)


class TestVcfExport(object):
    def test_protein_altering(self, vcf_writer, caplog, v600e):
        vcf_writer.addrecord(v600e)
        assert not caplog.records
        assert len(vcf_writer.variant_records) == 1
        out_dict = vcf_writer.writerecords()
        assert out_dict[0]['POS'] == '140453136'
        assert out_dict[0]['REF'] == 'A'
        assert out_dict[0]['ALT'] == 'T'

    def test_simple_insertion(self, vcf_writer, caplog, a56fs):
        assert a56fs.is_insertion
        vcf_writer.addrecord(a56fs)
        assert not caplog.records
        assert len(vcf_writer.variant_records) == 1
        out_dict = vcf_writer.writerecords()
        assert out_dict[0]['POS'] == '10183697'
        assert out_dict[0]['REF'] == 'G'
        assert out_dict[0]['ALT'] == 'GA'

    def test_simple_deletion(self, vcf_writer, caplog, v273fs):
        assert v273fs.is_deletion
        vcf_writer.addrecord(v273fs)
        assert not caplog.records
        assert len(vcf_writer.variant_records) == 1
        out_dict = vcf_writer.writerecords()
        assert out_dict[0]['POS'] == '47641432'
        assert out_dict[0]['REF'] == 'GT'
        assert out_dict[0]['ALT'] == 'G'

    def test_complex_insertion(self, vcf_writer, caplog, v2444fs):
        assert v2444fs.is_insertion
        vcf_writer.addrecord(v2444fs)
        assert not caplog.records
        assert len(vcf_writer.variant_records) == 1
        out_dict = vcf_writer.writerecords()
        assert out_dict[0]['POS'] == '139390861'
        assert out_dict[0]['REF'] == 'GG'
        assert out_dict[0]['ALT'] == 'GTGT'

    def test_complex_deletion(self, vcf_writer, caplog, l158fs):
        assert l158fs.is_deletion
        vcf_writer.addrecord(l158fs)
        assert not caplog.records
        assert len(vcf_writer.variant_records) == 1
        out_dict = vcf_writer.writerecords()
        assert out_dict[0]['POS'] == '10191480'
        assert out_dict[0]['REF'] == 'TGAA'
        assert out_dict[0]['ALT'] == 'TC'

    def test_addrecord_from_gene(self, vcf_writer):
        gene = civic.get_gene_by_id(24)
        vcf_writer.addrecord(gene)
        assert len(vcf_writer.variant_records) <= len(gene.variants)

    def test_addrecord_from_evidence(self, vcf_writer):
        evidence = civic._get_element_by_id('evidence', 373)
        vcf_writer.addrecord(evidence)
        assert len(vcf_writer.variant_records) == 1
        assert list(vcf_writer.variant_records)[0].id == evidence.variant.id

    def test_addrecords(self, vcf_writer, v600e, l158fs):
        vcf_writer.addrecords([v600e, l158fs])
        assert len(vcf_writer.variant_records) == 2

    def test_addrecord_wrong_type(self, vcf_writer):
        evidence = civic._get_element_by_id('evidence', 373)
        with pytest.raises(ValueError) as context:
            vcf_writer.addrecord(evidence.source)
        assert "Expected a CIViC Gene, Variant, Assertion or Evidence record" in str(context.value)
