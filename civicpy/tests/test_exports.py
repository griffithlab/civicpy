import pytest
from civicpy import exports
import io


@pytest.fixture(scope='module')
def vcf_stream():
    return io.StringIO()


@pytest.fixture(scope='module')
def vcf_writer(vcf_stream):
    return exports.VCFWriter(vcf_stream)


class TestVcfExport(object):

    def test_protein_altering(self, vcf_writer, caplog):
        assert False
        #TODO: get v600e from cache
        vcf_writer.addrecord(v600e)
        assert not caplog.records
        state1 = len(vcf_writer.evidence_records)
        assert state1 == len(v600e.evidence)
        vcf_writer.addrecord()
