import pytest
from civicpy import civic
from civicpy.exports.civic_vcf_record import CivicVcfRecord
import io


#snv
@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)

#simple insertion
@pytest.fixture(scope="module")
def a56fs():
    return civic.get_variant_by_id(1785)

#simple deletion
@pytest.fixture(scope="module")
def v273fs():
    return civic.get_variant_by_id(762)

#complex insertion
@pytest.fixture(scope="module")
def v2444fs():
    return civic.get_variant_by_id(137)

#complex deletion
@pytest.fixture(scope="module")
def l158fs():
    return civic.get_variant_by_id(2137)

class TestCivicVcfRecord(object):
    def test_protein_altering(self, caplog, v600e):
        record = CivicVcfRecord(v600e)
        assert not caplog.records
        assert record.POS == 140453136
        assert record.REF == 'A'
        assert record.ALT[0].value == 'T'

    def test_simple_insertion(self, caplog, a56fs):
        assert a56fs.is_insertion
        record = CivicVcfRecord(a56fs)
        assert not caplog.records
        assert record.POS == 10183697
        assert record.REF == 'G'
        assert record.ALT[0].value == 'GA'

    def test_simple_deletion(self, caplog, v273fs):
        assert v273fs.is_deletion
        record = CivicVcfRecord(v273fs)
        assert not caplog.records
        assert record.POS == 47641432
        assert record.REF == 'GT'
        assert record.ALT[0].value == 'G'

    def test_complex_insertion(self, caplog, v2444fs):
        assert v2444fs.is_insertion
        record = CivicVcfRecord(v2444fs)
        assert not caplog.records
        assert record.POS == 139390861
        assert record.REF == 'GG'
        assert record.ALT[0].value == 'GTGT'

    def test_complex_deletion(self, caplog, l158fs):
        assert l158fs.is_deletion
        record = CivicVcfRecord(l158fs)
        assert not caplog.records
        assert record.POS == 10191480
        assert record.REF == 'TGAA'
        assert record.ALT[0].value == 'TC'

    def test_addrecord_from_gene(self):
        gene = civic.get_gene_by_id(24)
        records = [CivicVcfRecord(v) for v in gene.variants if v.is_valid_for_vcf()]
        assert len(records) <= len(gene.variants)

    def test_addrecord_from_evidence(self):
        evidence = civic._get_element_by_id('evidence', 12)
        records = [CivicVcfRecord(v) for v in evidence.variants if v.is_valid_for_vcf()]
        assert len(records) == 1
        assert int(records[0].ID[0]) == evidence.molecular_profile.variants[0].id

    def test_addrecord_from_assertion(self):
        assertion = civic._get_element_by_id('assertion', 7)
        records = [CivicVcfRecord(v) for v in assertion.variants if v.is_valid_for_vcf()]
        assert len(records) == 1
        assert int(records[0].ID[0]) == assertion.molecular_profile.variants[0].id

    def test_add_record_from_molecular_profile(self):
        mp = civic._get_element_by_id('molecular_profile', 12)
        records = [CivicVcfRecord(v) for v in mp.variants if v.is_valid_for_vcf()]
        assert len(records) == 1
        assert int(records[0].ID[0]) == mp.variants[0].id

    def test_addrecord_wrong_type(self):
        evidence = civic._get_element_by_id('evidence', 373)
        with pytest.raises(Exception) as context:
            CivicVcfRecord(evidence)
        assert "Variant is not a GeneVariant" in str(context.value)
