from unittest.mock import PropertyMock, patch
import pytest
from deepdiff import DeepDiff
from ga4gh.va_spec.base import Condition, ConditionSet, Statement, TherapyGroup
from ga4gh.va_spec.aac_2017 import (
    VariantTherapeuticResponseStudyStatement,
    VariantDiagnosticStudyStatement,
    VariantPrognosticStudyStatement,
)


from civicpy import civic
from civicpy.exports.civic_vcf_record import CivicVcfRecord
from civicpy.exports.civic_gks_record import (
    CivicGksEvidence,
    CivicGksMolecularProfile,
    CivicGksRecordError,
    CivicGksPredictiveAssertion,
    CivicGksDiagnosticAssertion,
    CivicGksPrognosticAssertion,
    CivicGksTherapyGroup,
    create_gks_record_from_assertion,
)


# snv
@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)


@pytest.fixture(scope="module")
def v600e_mp():
    return civic.get_molecular_profile_by_id(12)


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
def eid9285():
    """Create test fixture for functional evidence"""
    return civic.get_evidence_by_id(9285)


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
def aid93():
    """Create test fixture for assertion with single phenotype"""
    return civic.get_assertion_by_id(93)


@pytest.fixture(scope="module")
def aid115():
    """Create test fixture for assertion with phenotypes"""
    return civic.get_assertion_by_id(115)


@pytest.fixture(scope="module")
def aid117():
    """Create test fixture for assertion not supported for GKS"""
    return civic.get_assertion_by_id(117)


@pytest.fixture(scope="module")
def approval4():
    """Create test fixture for active approval"""
    return civic.get_approval_by_id(4)


@pytest.fixture(scope="module")
def gks_contributions():
    return [
        {
            "type": "Contribution",
            "contributor": {
                "id": "civic.organization:1",
                "type": "Agent",
                "name": "The McDonnell Genome Institute",
                "description": "The McDonnell Genome Institute (MGI) is a world leader in the fast-paced, constantly changing field of genomics. A truly unique institution, we are pushing the limits of academic research by creating, testing, and implementing new approaches to the study of biology with the goal of understanding human health and disease, as well as evolution and the biology of other organisms.",
                "extensions": [{"name": "is_approved_vcep", "value": False}],
            },
            "activityType": "approval.last_reviewed",
            "date": "2025-05-27",
        }
    ]


@pytest.fixture(scope="module")
def gks_method():
    """Create test fixture for method GKS representation."""
    return {
        "id": "civic.method:2019",
        "name": "CIViC Curation SOP (2019)",
        "reportedIn": {
            "id": "pmid:31779674",
            "name": "Danos et al., 2019, Genome Med.",
            "title": "Standard operating procedure for curation and clinical interpretation of variants in cancer",
            "doi": "10.1186/s13073-019-0687-x",
            "pmid": "31779674",
            "urls": [
                "https://doi.org/10.1186/s13073-019-0687-x",
                "https://pubmed.ncbi.nlm.nih.gov/31779674/",
            ],
            "aliases": ["CIViC curation SOP"],
            "type": "Document",
        },
        "methodType": "curation",
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
        "aliases": ["LEU858ARG", "L813R", "LEU813ARG"],
        "mappings": [
            {
                "coding": {
                    "code": "CA126713",
                    "system": "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "code": "16609",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "code": "376282",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "code": "376280",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "code": "rs121434568",
                    "system": "https://www.ncbi.nlm.nih.gov/snp/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "code": "33",
                    "id": "civic.mpid:33",
                    "system": "https://civicdb.org/links/molecular_profile/",
                },
                "relation": "exactMatch",
            },
            {
                "coding": {
                    "code": "33",
                    "id": "civic.vid:33",
                    "name": "L858R",
                    "system": "https://civicdb.org/links/variant/",
                    "extensions": [
                        {"name": "subtype", "value": "gene_variant"},
                        {
                            "name": "variant_types",
                            "value": [
                                {
                                    "coding": {
                                        "id": "civic.variant_type:47",
                                        "code": "SO:0001583",
                                        "name": "Missense Variant",
                                        "system": "http://www.sequenceontology.org/browser/current_svn/term/",
                                    },
                                    "relation": "exactMatch",
                                }
                            ],
                        },
                    ],
                },
                "relation": "exactMatch",
            },
        ],
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
        "id": "civic.sid:1725",
        "name": "Dungo et al., 2013",
        "title": "Afatinib: first global approval.",
        "pmid": "23982599",
        "type": "Document",
        "urls": [
            "https://civicdb.org/links/source/1725",
            "http://www.ncbi.nlm.nih.gov/pubmed/23982599",
        ],
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
            "mappings": [
                {
                    "coding": {
                        "code": "e000001",
                        "name": "authoritative evidence",
                        "system": "https://go.osu.edu/evidence-codes",
                    },
                    "relation": "exactMatch",
                }
            ],
        },
        "proposition": {
            "type": "VariantTherapeuticResponseProposition",
            "predicate": "predictsSensitivityTo",
            "objectTherapeutic": gks_tid146,
            "conditionQualifier": gks_did8,
            "alleleOriginQualifier": {
                "name": "somatic",
                "extensions": [{"name": "civic_variant_origin", "value": "SOMATIC"}],
            },
            "geneContextQualifier": gks_gid19,
            "subjectVariant": gks_mpid33,
        },
        "specifiedBy": gks_method,
        "reportedIn": [gks_source592, "https://civicdb.org/links/evidence/2997"],
    }
    return Statement(**params)


@pytest.fixture(scope="module")
def gks_aid6(
    gks_contributions,
    gks_method,
    gks_mpid33,
    gks_gid19,
    gks_tid146,
    gks_did8,
    gks_eid2997,
):
    """Create CIVIC AID6 GKS representation."""
    params = {
        "id": "civic.aid:6",
        "contributions": gks_contributions,
        "description": "L858R is among the most common sensitizing EGFR mutations in NSCLC, and is assessed via DNA mutational analysis, including Sanger sequencing and next generation sequencing methods. Tyrosine kinase inhibitor afatinib is FDA approved as a first line systemic therapy in NSCLC with sensitizing EGFR mutation (civic.EID:2997).",
        "type": "Statement",
        "specifiedBy": gks_method,
        "proposition": {
            "type": "VariantTherapeuticResponseProposition",
            "subjectVariant": gks_mpid33,
            "geneContextQualifier": gks_gid19,
            "alleleOriginQualifier": {
                "name": "somatic",
                "extensions": [{"name": "civic_variant_origin", "value": "SOMATIC"}],
            },
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
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "Level A",
                        "system": "AMP/ASCO/CAP (AAC) Guidelines, 2017",
                    },
                },
            }
        ],
        "reportedIn": ["https://civicdb.org/links/assertion/6"],
    }
    return VariantTherapeuticResponseStudyStatement(**params)


@pytest.fixture(scope="module")
def gks_aid93_object_condition():
    """Create test fixture for GKS AID 93 object condition"""
    return {
        "conditions": [
            {
                "id": "civic.did:3225",
                "conceptType": "Disease",
                "name": "CNS Neuroblastoma With FOXR2 Activation",
                "mappings": [
                    {
                        "coding": {
                            "system": "https://disease-ontology.org/?id=",
                            "code": "DOID:0080906",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
            {
                "id": "civic.phenotype:15320",
                "conceptType": "Phenotype",
                "name": "Pediatric onset",
                "mappings": [
                    {
                        "coding": {
                            "system": "https://hpo.jax.org/app/browse/term/",
                            "code": "HP:0410280",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        ],
        "membershipOperator": "AND",
    }


@pytest.fixture(scope="module")
def gks_aid115_object_condition():
    """Create test fixture for GKS AID 115 object condition"""
    return {
        "conditions": [
            {
                "id": "civic.did:3387",
                "conceptType": "Disease",
                "name": "Diffuse Astrocytoma, MYB- Or MYBL1-altered",
                "mappings": [
                    {
                        "coding": {
                            "system": "https://disease-ontology.org/?id=",
                            "code": "DOID:0081279",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
            {
                "conditions": [
                    {
                        "id": "civic.phenotype:8121",
                        "conceptType": "Phenotype",
                        "name": "Childhood onset",
                        "mappings": [
                            {
                                "coding": {
                                    "system": "https://hpo.jax.org/app/browse/term/",
                                    "code": "HP:0011463",
                                },
                                "relation": "exactMatch",
                            }
                        ],
                    },
                    {
                        "id": "civic.phenotype:2656",
                        "conceptType": "Phenotype",
                        "name": "Juvenile onset",
                        "mappings": [
                            {
                                "coding": {
                                    "system": "https://hpo.jax.org/app/browse/term/",
                                    "code": "HP:0003621",
                                },
                                "relation": "exactMatch",
                            }
                        ],
                    },
                    {
                        "id": "civic.phenotype:2643",
                        "conceptType": "Phenotype",
                        "name": "Adult onset",
                        "mappings": [
                            {
                                "coding": {
                                    "system": "https://hpo.jax.org/app/browse/term/",
                                    "code": "HP:0003581",
                                },
                                "relation": "exactMatch",
                            }
                        ],
                    },
                ],
                "membershipOperator": "OR",
            },
        ],
        "membershipOperator": "AND",
    }


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


class TestCivicGksMolecularProfile(object):
    """Test that CivicGksMolecularProfile works as expected"""

    def test_get_extensions(self, v600e_mp):
        """Test that get_extensions method works as expected"""
        variant = v600e_mp.variants[0]

        with patch.object(
            variant, "hgvs_expressions", new=["N/A", "XR_001744858.1:n.1823-3918T>A"]
        ):
            gks_mp = CivicGksMolecularProfile(v600e_mp)
            extensions = gks_mp.get_extensions(v600e_mp)
            assert extensions

            expressions = next(
                (ext for ext in extensions if ext.name == "expressions"), None
            )
            assert expressions is None

    def test_na_clinvar_mapping(self, v600e_mp):
        """Test that get_aliases_and_mappings method works as expected when no clinvar entry found"""
        variant = v600e_mp.variants[0]

        with patch.object(variant, "clinvar_entries", new=["N/A"]):
            gks_mp = CivicGksMolecularProfile(v600e_mp)
            _, mappings = gks_mp.get_aliases_and_mappings(v600e_mp)
            assert mappings
            assert not any(
                m.coding.system == "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
                for m in mappings
            )


class TestCivicGksTherapyGroup(object):
    """Test that CivicGksTherapyGroup works as expected"""

    def test_no_therapies(self):
        """Test that CivicGksTherapyGroup works as expected when no therapies provided"""
        with pytest.raises(CivicGksRecordError, match=r"No therapies provided"):
            CivicGksTherapyGroup(therapies=[], therapy_interaction_type=None)


class TestCivicGksEvidence(object):
    """Test that CivicGksEvidence works as expected"""

    def test_invalid(self, eid9285):
        """Test that invalid assertions raises custom exception"""
        with pytest.raises(
            CivicGksRecordError, match=r"Evidence 9285 is not valid for GKS."
        ):
            CivicGksEvidence(eid9285)


class TestCivicGksPredictiveAssertion(object):
    """Test that CivicGksPredictiveAssertion works as expected"""

    def test_valid_single_therapy(self, aid6, approval4, gks_aid6):
        """Test that single therapy works as expected"""
        record = CivicGksPredictiveAssertion(aid6, approval=approval4)
        assert isinstance(record, VariantTherapeuticResponseStudyStatement)
        assert len(record.hasEvidenceLines) == 1

        # Don't need to test ALL has evidence lines
        check_evs = []
        el = record.hasEvidenceLines[0]
        assert len(el.hasEvidenceItems) == 6
        for ev in el.hasEvidenceItems:
            if ev.id == "civic.eid:2997":
                check_evs.append(ev)

        record_copy = record.model_copy(deep=True)
        record_copy.hasEvidenceLines[0].hasEvidenceItems = check_evs
        record_copy = record_copy.model_dump(exclude_none=True)
        diff = DeepDiff(
            record_copy, gks_aid6.model_dump(exclude_none=True), ignore_order=True
        )
        assert diff == {}, gks_aid6.id

    def test_valid_combination_therapy(self, aid7):
        """Test that combination therapy works as expected"""
        record = CivicGksPredictiveAssertion(aid7)
        assert isinstance(record, VariantTherapeuticResponseStudyStatement)
        assert len(record.hasEvidenceLines) == 1
        assert len(record.hasEvidenceLines[0].hasEvidenceItems) == 4
        therapy = record.proposition.objectTherapeutic.root
        assert isinstance(therapy, TherapyGroup)
        assert therapy.membershipOperator == "AND"
        assert len(therapy.therapies) == 2
        therapy_ids = {t.id for t in therapy.therapies}
        assert therapy_ids == {"civic.tid:19", "civic.tid:22"}

    @patch.object(civic.Assertion, "is_valid_for_gks_json")
    @patch.object(civic.Assertion, "evidence_items")
    @patch.object(civic.FusionVariant, "hgvs_expressions", create=True)
    @patch.object(
        civic.FusionVariant,
        "allele_registry_id",
        create=True,
        new_callable=PropertyMock,
    )
    @patch.object(
        civic.FusionVariant, "clinvar_entries", create=True, new_callable=PropertyMock
    )
    @patch.object(civic.FusionVariant, "coordinates", create=True)
    @patch.object(
        civic.FusionVariant, "gene", new=civic.get_gene_by_id(1590), create=True
    )
    def test_valid_substitution_therapy(
        self,
        test_coordinates,
        test_clinvar_entries,
        test_allele_registry_id,
        test_hgvs_expressions,
        test_evidence_items,
        test_is_valid_for_gks_json,
        aid19,
    ):
        """Test that substitution therapy works as expected"""
        test_coordinates.return_value = None
        test_clinvar_entries.return_value = []
        test_allele_registry_id.return_value = None
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
        assert len(record.hasEvidenceLines) == 1
        assert len(record.hasEvidenceLines[0].hasEvidenceItems) == 6
        assert record.proposition.predicate == "associatedWithWorseOutcomeFor"
        assert record.strength.primaryCoding.code.root == "Level A"
        assert record.classification.primaryCoding.code.root == "Tier I"


class TestCivicGksDiagnosticAssertion(object):
    """Test that CivicGksDiagnosticAssertion works as expected"""

    def test_valid(
        self,
        aid9,
        aid93,
        gks_aid93_object_condition,
        aid115,
        gks_aid115_object_condition,
    ):
        """Test that valid assertion works as expected"""
        record = CivicGksDiagnosticAssertion(aid9)
        assert isinstance(record, VariantDiagnosticStudyStatement)
        assert len(record.hasEvidenceLines) == 1
        assert len(record.hasEvidenceLines[0].hasEvidenceItems) == 2
        assert record.proposition.predicate == "isDiagnosticInclusionCriterionFor"
        assert record.strength.primaryCoding.code.root == "Level C"
        assert record.classification.primaryCoding.code.root == "Tier II"

        # Single phenotype (complex condition set)
        record = CivicGksDiagnosticAssertion(aid93)
        assert isinstance(record, VariantDiagnosticStudyStatement)
        record_object_condition = record.proposition.objectCondition
        assert isinstance(record_object_condition, Condition)
        assert isinstance(record_object_condition.root, ConditionSet)
        diff = DeepDiff(
            record_object_condition.model_dump(exclude_none=True),
            gks_aid93_object_condition,
            ignore_order=True,
        )
        assert diff == {}

        # Phenotypes (complex condition set)
        record = CivicGksDiagnosticAssertion(aid115)
        assert isinstance(record, VariantDiagnosticStudyStatement)
        record_object_condition = record.proposition.objectCondition
        assert isinstance(record_object_condition, Condition)
        assert isinstance(record_object_condition.root, ConditionSet)
        diff = DeepDiff(
            record_object_condition.model_dump(exclude_none=True),
            gks_aid115_object_condition,
            ignore_order=True,
        )
        assert diff == {}


class TestCivicGksRecord(object):
    """Test that GKS Record helper methods work correctly"""

    def test_invalid(self, aid117):
        """Test that unsupported assertion types raise NotImplementedError"""

        with pytest.raises(
            NotImplementedError,
            match=r"Assertion type ONCOGENIC is not currently supported",
        ):
            create_gks_record_from_assertion(aid117)
