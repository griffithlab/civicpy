from unittest.mock import patch
import pytest
from deepdiff import DeepDiff
from ga4gh.va_spec.base import Statement, EvidenceLine, TherapyGroup
from ga4gh.va_spec.aac_2017 import (
    VariantTherapeuticResponseStudyStatement,
    VariantDiagnosticStudyStatement,
    VariantPrognosticStudyStatement,
)


from civicpy import civic
from civicpy.exports.civic_vcf_record import CivicVcfRecord
from civicpy.exports.civic_gks_record import (
    CivicGksRecordError,
    CivicGksPredictiveAssertion,
    CivicGksDiagnosticAssertion,
    CivicGksPrognosticAssertion,
)


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


@pytest.fixture(scope="module")
def aid6():
    """Create test fixture for predictive assertion (single therapy)"""
    return civic.get_assertion_by_id(6)


@pytest.fixture(scope="module")
def aid7():
    """Create test fixture for predictive assertion (combination therapy)"""
    return civic.get_assertion_by_id(7)


@pytest.fixture(scope="module")
def aid9():
    """Create test fixture for diagnostic assertion"""
    return civic.get_assertion_by_id(9)


@pytest.fixture(scope="module")
def aid19():
    """Create test fixture for predictive assertion (substitution therapy)"""
    return civic.get_assertion_by_id(19)


@pytest.fixture(scope="module")
def aid20():
    """Create test fixture for prognostic assertion"""
    return civic.get_assertion_by_id(20)


@pytest.fixture(scope="module")
def aid117():
    """Create test fixture for assertion not supported for GKS"""
    return civic.get_assertion_by_id(117)


@pytest.fixture(scope="module")
def gks_method():
    """Create test fixture for method GKS representation."""
    return {
        "id": "civic.method:2019",
        "name": "CIViC Curation SOP (2019)",
        "reportedIn": {
            "name": "Danos et al., 2019, Genome Med.",
            "title": "Standard operating procedure for curation and clinical interpretation of variants in cancer",
            "doi": "10.1186/s13073-019-0687-x",
            "pmid": 31779674,
            "type": "Document",
        },
        "methodType": "variant curation standard operating procedure",
        "type": "Method",
    }


@pytest.fixture(scope="module")
def gks_mpid33():
    """Create CIViC MPID 33 GKS representation"""
    return {
        "id": "civic.mpid:33",
        "type": "CategoricalVariant",
        "description": "EGFR L858R has long been recognized as a functionally significant mutation in cancer, and is one of the most prevalent single mutations in lung cancer. Best described in non-small cell lung cancer (NSCLC), the mutation seems to confer sensitivity to first and second generation TKI's like gefitinib and neratinib. NSCLC patients with this mutation treated with TKI's show increased overall and progression-free survival, as compared to chemotherapy alone. Third generation TKI's are currently in clinical trials that specifically focus on mutant forms of EGFR, a few of which have shown efficacy in treating patients that failed to respond to earlier generation TKI therapies.",
        "name": "EGFR L858R",
        "aliases": ["LEU858ARG", "L813R", "LEU813ARG", "RS121434568"],
        "extensions": [
            {
                "name": "CIViC representative coordinate",
                "value": {
                    "chromosome": "7",
                    "start": 55259515,
                    "stop": 55259515,
                    "reference_bases": "T",
                    "variant_bases": "G",
                    "representative_transcript": "ENST00000275493.2",
                    "ensembl_version": 75,
                    "reference_build": "GRCh37",
                    "type": "coordinates",
                },
            },
            {
                "name": "CIViC Molecular Profile Score",
                "value": 379.0,
            },
            {
                "name": "expressions",
                "value": [
                    {"syntax": "hgvs.c", "value": "ENST00000275493.2:c.2573T>G"},
                    {"syntax": "hgvs.c", "value": "NM_005228.4:c.2573T>G"},
                    {"syntax": "hgvs.g", "value": "NC_000007.13:g.55259515T>G"},
                    {"syntax": "hgvs.p", "value": "NP_005219.2:p.Leu858Arg"},
                ],
            },
        ],
    }


@pytest.fixture(scope="module")
def gks_gid19():
    """Create test fixture for CIViC GID19 GKS representation."""
    return {
        "id": "civic.gid:19",
        "conceptType": "Gene",
        "name": "EGFR",
        "mappings": [
            {
                "coding": {
                    "id": "ncbigene:1956",
                    "code": "1956",
                    "system": "https://www.ncbi.nlm.nih.gov/gene/",
                },
                "relation": "exactMatch",
            },
        ],
        "extensions": [
            {
                "name": "description",
                "value": "EGFR is widely recognized for its importance in cancer. Amplification and mutations have been shown to be driving events in many cancer types. Its role in non-small cell lung cancer, glioblastoma and basal-like breast cancers has spurred many research and drug development efforts. Tyrosine kinase inhibitors have shown efficacy in EGFR amplfied tumors, most notably gefitinib and erlotinib. Mutations in EGFR have been shown to confer resistance to these drugs, particularly the variant T790M, which has been functionally characterized as a resistance marker for both of these drugs. The later generation TKI's have seen some success in treating these resistant cases, and targeted sequencing of the EGFR locus has become a common practice in treatment of non-small cell lung cancer. Overproduction of ligands is another possible mechanism of activation of EGFR. ERBB ligands include EGF, TGF-a, AREG, EPG, BTC, HB-EGF, EPR and NRG1-4 (for detailed information please refer to the respective ligand section).",
            },
            {
                "name": "aliases",
                "value": [
                    "EGFR",
                    "ERBB",
                    "ERBB1",
                    "ERRP",
                    "HER1",
                    "NISBD2",
                    "PIG61",
                    "mENA",
                ],
            },
        ],
    }


@pytest.fixture(scope="module")
def gks_did8():
    """Create test fixture for CIViC DID8 GKS representation."""
    return {
        "id": "civic.did:8",
        "conceptType": "Disease",
        "name": "Lung Non-small Cell Carcinoma",
        "mappings": [
            {
                "coding": {
                    "id": "DOID:3908",
                    "code": "DOID:3908",
                    "system": "https://disease-ontology.org/?id=",
                },
                "relation": "exactMatch",
            },
        ],
    }


@pytest.fixture(scope="module")
def gks_tid146():
    """Create test fixture for CIViC TID146 GKS representation."""
    return {
        "id": "civic.tid:146",
        "conceptType": "Therapy",
        "name": "Afatinib",
        "mappings": [
            {
                "coding": {
                    "id": "ncit:C66940",
                    "code": "C66940",
                    "system": "https://ncit.nci.nih.gov/ncitbrowser/ConceptReport.jsp?dictionary=NCI_Thesaurus&code=",
                },
                "relation": "exactMatch",
            },
        ],
        "extensions": [
            {
                "name": "aliases",
                "value": [
                    "(2e)-N-(4-(3-Chloro-4-Fluoroanilino)-7-(((3s)-Oxolan-3-yl)Oxy)Quinoxazolin-6-yl)-4-(Dimethylamino)But-2-Enamide",
                    "BIBW 2992",
                    "BIBW-2992",
                    "BIBW2992",
                ],
            },
        ],
    }


@pytest.fixture(scope="module")
def gks_source592():
    """Create fixture for source 592 GKS representation"""
    return {
        "id": "civic.source:1725",
        "name": "Dungo et al., 2013",
        "title": "Afatinib: first global approval.",
        "pmid": 23982599,
        "type": "Document",
    }


@pytest.fixture(scope="module")
def gks_eid2997(
    gks_mpid33,
    gks_tid146,
    gks_did8,
    gks_gid19,
    gks_method,
    gks_source592,
):
    """Create CIVIC EID2997 GKS representation."""
    params = {
        "id": "civic.eid:2997",
        "type": "Statement",
        "description": "Afatinib, an irreversible inhibitor of the ErbB family of tyrosine kinases has been approved in the US for the first-line treatment of patients with metastatic non-small-cell lung cancer (NSCLC) who have tumours with EGFR exon 19 deletions or exon 21 (L858R) substitution mutations as detected by a US FDA-approved test",
        "direction": "supports",
        "strength": {
            "name": "Validated association",
            "primaryCoding": {
                "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                "code": "A",
            },
        },
        "proposition": {
            "type": "VariantTherapeuticResponseProposition",
            "predicate": "predictsSensitivityTo",
            "objectTherapeutic": gks_tid146,
            "conditionQualifier": gks_did8,
            "alleleOriginQualifier": {"name": "SOMATIC"},
            "geneContextQualifier": gks_gid19,
            "subjectVariant": gks_mpid33,
        },
        "specifiedBy": gks_method,
        "reportedIn": [gks_source592],
    }
    return Statement(**params)


@pytest.fixture(scope="module")
def gks_aid6(gks_method, gks_mpid33, gks_gid19, gks_tid146, gks_did8, gks_eid2997):
    """Create CIVIC AID6 GKS representation."""
    params = {
        "id": "civic.aid:6",
        "description": "L858R is among the most common sensitizing EGFR mutations in NSCLC, and is assessed via DNA mutational analysis, including Sanger sequencing and next generation sequencing methods. Tyrosine kinase inhibitor afatinib is FDA approved as a first line systemic therapy in NSCLC with sensitizing EGFR mutation (civic.EID:2997).",
        "type": "Statement",
        "specifiedBy": gks_method,
        "proposition": {
            "type": "VariantTherapeuticResponseProposition",
            "subjectVariant": gks_mpid33,
            "geneContextQualifier": gks_gid19,
            "alleleOriginQualifier": {"name": "SOMATIC"},
            "predicate": "predictsSensitivityTo",
            "objectTherapeutic": gks_tid146,
            "conditionQualifier": gks_did8,
        },
        "direction": "supports",
        "strength": {
            "primaryCoding": {
                "system": "AMP/ASCO/CAP (AAC) Guidelines, 2017",
                "code": "Level A",
            },
            "mappings": [
                {
                    "coding": {
                        "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                        "code": "A",
                        "name": "Validated association",
                    },
                    "relation": "exactMatch",
                }
            ],
        },
        "classification": {
            "primaryCoding": {
                "system": "AMP/ASCO/CAP (AAC) Guidelines, 2017",
                "code": "Tier I",
            },
        },
        "hasEvidenceLines": [
            {
                "type": "EvidenceLine",
                "hasEvidenceItems": [gks_eid2997],
                "directionOfEvidenceProvided": "supports",
            }
        ],
    }
    return VariantTherapeuticResponseStudyStatement(**params)


class TestCivicVcfRecord(object):
    def test_protein_altering(self, caplog, v600e):
        record = CivicVcfRecord(v600e)
        assert not caplog.records
        assert record.POS == 140453136
        assert record.REF == "A"
        assert record.ALT[0].value == "T"

    def test_simple_insertion(self, caplog, a56fs):
        assert a56fs.is_insertion
        record = CivicVcfRecord(a56fs)
        assert not caplog.records
        assert record.POS == 10183697
        assert record.REF == "G"
        assert record.ALT[0].value == "GA"

    def test_simple_deletion(self, caplog, v273fs):
        assert v273fs.is_deletion
        record = CivicVcfRecord(v273fs)
        assert not caplog.records
        assert record.POS == 47641432
        assert record.REF == "GT"
        assert record.ALT[0].value == "G"

    def test_complex_insertion(self, caplog, v2444fs):
        assert v2444fs.is_insertion
        record = CivicVcfRecord(v2444fs)
        assert not caplog.records
        assert record.POS == 139390861
        assert record.REF == "GG"
        assert record.ALT[0].value == "GTGT"

    def test_complex_deletion(self, caplog, l158fs):
        assert l158fs.is_deletion
        record = CivicVcfRecord(l158fs)
        assert not caplog.records
        assert record.POS == 10191480
        assert record.REF == "TGAA"
        assert record.ALT[0].value == "TC"

    def test_addrecord_from_gene(self):
        gene = civic.get_gene_by_id(24)
        records = [CivicVcfRecord(v) for v in gene.variants if v.is_valid_for_vcf()]
        assert len(records) <= len(gene.variants)

    def test_addrecord_from_evidence(self):
        evidence = civic._get_element_by_id("evidence", 12)
        records = [CivicVcfRecord(v) for v in evidence.variants if v.is_valid_for_vcf()]
        assert len(records) == 1
        assert int(records[0].ID[0]) == evidence.molecular_profile.variants[0].id

    def test_addrecord_from_assertion(self):
        assertion = civic._get_element_by_id("assertion", 7)
        records = [
            CivicVcfRecord(v) for v in assertion.variants if v.is_valid_for_vcf()
        ]
        assert len(records) == 1
        assert int(records[0].ID[0]) == assertion.molecular_profile.variants[0].id

    def test_add_record_from_molecular_profile(self):
        mp = civic._get_element_by_id("molecular_profile", 12)
        records = [CivicVcfRecord(v) for v in mp.variants if v.is_valid_for_vcf()]
        assert len(records) == 1
        assert int(records[0].ID[0]) == mp.variants[0].id

    def test_addrecord_wrong_type(self):
        evidence = civic._get_element_by_id("evidence", 373)
        with pytest.raises(Exception) as context:
            CivicVcfRecord(evidence)
        assert "Variant is not a GeneVariant" in str(context.value)


class TestCivicGksPredictiveAssertion(object):
    """Test that CivicGksPredictiveAssertion works as expected"""

    def test_valid_single_therapy(self, aid6, gks_aid6):
        """Test that single therapy works as expected"""
        record = CivicGksPredictiveAssertion(aid6)
        assert isinstance(record, VariantTherapeuticResponseStudyStatement)
        assert len(record.hasEvidenceLines) > 1

        # Don't need to test ALL has evidence lines
        record_copy = record.model_copy(deep=True)
        check_els = []
        for el in record_copy.hasEvidenceLines:
            assert isinstance(el, EvidenceLine)
            assert len(el.hasEvidenceItems) == 1

            if el.hasEvidenceItems[0].id == "civic.eid:2997":
                check_els.append(el)

        record_copy.hasEvidenceLines = check_els
        record_copy = record_copy.model_dump(exclude_none=True)
        diff = DeepDiff(
            record_copy, gks_aid6.model_dump(exclude_none=True), ignore_order=True
        )
        assert diff == {}, gks_aid6.id

    def test_valid_combination_therapy(self, aid7):
        """Test that combination therapy works as expected"""
        record = CivicGksPredictiveAssertion(aid7)
        assert isinstance(record, VariantTherapeuticResponseStudyStatement)
        assert len(record.hasEvidenceLines) > 1
        therapy = record.proposition.objectTherapeutic.root
        assert isinstance(therapy, TherapyGroup)
        assert therapy.membershipOperator == "AND"
        assert len(therapy.therapies) == 2
        therapy_ids = {t.id for t in therapy.therapies}
        assert therapy_ids == {"civic.tid:19", "civic.tid:22"}

    @patch.object(civic.Assertion, "is_valid_for_gks_json")
    @patch.object(civic.Assertion, "evidence_items")
    @patch.object(civic.FusionVariant, "hgvs_expressions", create=True)
    @patch.object(civic.FusionVariant, "coordinates", create=True)
    @patch.object(
        civic.FusionVariant, "gene", new=civic.get_gene_by_id(1590), create=True
    )
    def test_valid_substitution_therapy(
        self,
        test_coordinates,
        test_hgvs_expressions,
        test_evidence_items,
        test_is_valid_for_gks_json,
        aid19,
    ):
        """Test that substitution therapy works as expected"""
        test_coordinates.return_value = None
        test_is_valid_for_gks_json.return_value = True
        test_evidence_items.return_value = []
        test_hgvs_expressions.return_value = None
        record = CivicGksPredictiveAssertion(aid19)
        assert isinstance(record, VariantTherapeuticResponseStudyStatement)
        therapy = record.proposition.objectTherapeutic.root
        assert isinstance(therapy, TherapyGroup)
        assert therapy.membershipOperator == "OR"
        assert len(therapy.therapies) == 2
        therapy_ids = {t.id for t in therapy.therapies}
        assert therapy_ids == {"civic.tid:5", "civic.tid:20"}

    def test_invalid(self, aid117):
        """Test that invalid assertions raises custom exception"""
        with pytest.raises(
            CivicGksRecordError, match=r"Assertion is not valid for GKS."
        ):
            CivicGksPredictiveAssertion(aid117)


class TestCivicGksPrognosticAssertion(object):
    """Test that CivicGksPrognosticAssertion works as expected"""

    def test_valid(self, aid20):
        """Test that valid assertion works as expected"""
        record = CivicGksPrognosticAssertion(aid20)
        assert isinstance(record, VariantPrognosticStudyStatement)
        assert len(record.hasEvidenceLines) > 1
        assert record.proposition.predicate == "associatedWithWorseOutcomeFor"
        assert record.strength.primaryCoding.code.root == "Level A"
        assert record.classification.primaryCoding.code.root == "Tier I"


class TestCivicGksDiagnosticAssertion(object):
    """Test that CivicGksDiagnosticAssertion works as expected"""

    def test_valid(self, aid9):
        """Test that valid assertion works as expected"""
        record = CivicGksDiagnosticAssertion(aid9)
        assert isinstance(record, VariantDiagnosticStudyStatement)
        assert len(record.hasEvidenceLines) > 1
        assert record.proposition.predicate == "isDiagnosticInclusionCriterionFor"
        assert record.strength.primaryCoding.code.root == "Level C"
        assert record.classification.primaryCoding.code.root == "Tier II"
