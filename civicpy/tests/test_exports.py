import pytest
from civicpy import exports, civic
import io


@pytest.fixture(scope='module')
def vcf_stream():
    return io.StringIO()


@pytest.fixture(scope='module')
def vcf_writer(vcf_stream):
    return exports.VCFWriter(vcf_stream)

@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)

class TestVcfExport(object):

    #@pytest.mark.skip(reason="Implementation under development")
    def test_protein_altering(self, vcf_writer, caplog, v600e):
        vcf_writer.addrecord(v600e)
        assert not caplog.records
        assert len(vcf_writer.variant_records) == 1
        vcf_writer.writerecords()
