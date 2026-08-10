import re
from copy import deepcopy
from unittest.mock import Mock, PropertyMock, patch

import pytest
from deepdiff import DeepDiff
from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import Extension, iriReference
from ga4gh.va_spec.aac_2017 import (
    VariantClinicalSignificanceStatement,
)
from ga4gh.va_spec.base import Condition, ConditionSet, Statement, TherapyGroup
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
from ga4gh.vrs.models import Allele, CopyNumberCount, Expression, Syntax

from civicpy import civic
from civicpy.exports.civic_gks_record import (
    CivicGksClinSigAssertion,
    CivicGksEvidence,
    CivicGksGene,
    CivicGksMolecularProfile,
    CivicGksOncogenicAssertion,
    CivicGksRecordError,
    CivicGksTherapyGroup,
    ClinVarSubmissionType,
    create_gks_record_from_assertion,
)
from civicpy.exports.civic_vcf_record import CivicVcfRecord
from civicpy.exports.variation_normalizer import VariationNormalizerDataProxy


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
def aid202():
    """Create test fixture for oncogenic assertion"""
    return civic.get_assertion_by_id(202)


@pytest.fixture(scope="module")
def gks_contributions():
    return [
        {
            "type": "Contribution",
            "contributor": {
                "id": "civic.organization:1",
                "type": "Agent",
                "name": "CIViC",
                "description": "The CIViC Organization (formerly “The McDonnell Genome Institute” CIViC organization) comprises the founders, developers, editors, curators, and administrators who build and maintain the knowledgebase, based at Washington University in St. Louis. This group is dedicated to ensuring that high-quality cancer variant interpretations are broadly accessible for precision oncology. One of their main roles is evaluating and synthesizing crowdsourced community contributions into formal clinical Assertions. Once these Assertions meet the strict criteria of the CIViC standard operating procedure, core approval members approve them for 1-star submission to the ClinVar CIViC organization.",
                "extensions": [{"name": "is_approved_vcep", "value": False}],
            },
            "activityType": "approval.last_reviewed",
            "date": "2026-04-24",
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
        "constraints": [
            {
                "type": "DefiningAlleleConstraint",
                "allele": {
                    "id": "ga4gh:VA.S41CcMJT2bcd8R4-qXZWH1PoHWNtG2PZ",
                    "type": "Allele",
                    "name": "EGFR L858R",
                    "digest": "S41CcMJT2bcd8R4-qXZWH1PoHWNtG2PZ",
                    "location": {
                        "id": "ga4gh:SL.v0_edynH98OIu-0QPVT5anCSOriAFSDQ",
                        "type": "SequenceLocation",
                        "digest": "v0_edynH98OIu-0QPVT5anCSOriAFSDQ",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                        },
                        "start": 857,
                        "end": 858,
                        "sequence": "L",
                    },
                    "state": {"type": "LiteralSequenceExpression", "sequence": "R"},
                    "expressions": [
                        {"syntax": "hgvs.p", "value": "NP_005219.2:p.Leu858Arg"}
                    ],
                },
                "relations": [
                    {
                        "primaryCoding": {
                            "system": "ga4gh-gks-term:allele-relation",
                            "code": "liftover_to",
                        }
                    },
                    {
                        "primaryCoding": {
                            "system": "http://www.sequenceontology.org",
                            "code": "translation_of",
                        }
                    },
                ],
            }
        ],
        "aliases": ["LEU858ARG", "L813R", "LEU813ARG"],
        "mappings": [
            {
                "coding": {
                    "id": "clingen.allele:CA126713",
                    "code": "CA126713",
                    "system": "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "clinvar:16609",
                    "code": "16609",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "clinvar:376282",
                    "code": "376282",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "clinvar:376280",
                    "code": "376280",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "dbsnp:rs121434568",
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
                "name": "hgvsDescriptions",
                "value": [
                    "ENST00000275493.2:c.2573T>G",
                    "NM_005228.4:c.2573T>G",
                    "NC_000007.13:g.55259515T>G",
                    "NC_000007.14:g.55191822T>G",
                    "NP_005219.2:p.Leu858Arg",
                    "ENSP00000275493.2:p.Leu858Arg",
                ],
            },
            {"name": "maneSelectTranscript", "value": "ENST00000275493.7:c.2573T>G"},
            {"name": "categoricalVariationType", "value": "ProteinSequenceConsequence"},
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
                    "NNCIS",
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
def gks_therapeutic_proposition(gks_mpid33, gks_gid19, gks_tid146, gks_did8):
    """Create test fixture for GKS therapeutic proposition"""
    return {
        "type": "VariantTherapeuticResponseProposition",
        "subjectVariant": gks_mpid33,
        "geneContextQualifier": gks_gid19,
        "alleleOriginQualifier": {
            "name": "somatic",
            "mappings": [
                {
                    "coding": {
                        "code": "SOMATIC",
                        "system": "https://civicdb.org",
                        "iris": [
                            "https://civic.readthedocs.io/en/latest/model/evidence/origin.html"
                        ],
                    },
                    "relation": "exactMatch",
                }
            ],
        },
        "predicate": "predictsSensitivityTo",
        "objectTherapeutic": gks_tid146,
        "conditionQualifier": gks_did8,
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
            "https://civicdb.org/links/evidence/2997",
            "https://civicdb.org/links/source/1725",
            "http://www.ncbi.nlm.nih.gov/pubmed/23982599",
        ],
    }


@pytest.fixture(scope="module")
def gks_eid2997(
    gks_therapeutic_proposition,
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
        "proposition": gks_therapeutic_proposition,
        "specifiedBy": gks_method,
        "reportedIn": [gks_source592],
    }
    return Statement(**params)


@pytest.fixture(scope="module")
def gks_aid6(gks_method, gks_therapeutic_proposition, gks_eid2997, gks_source592):
    """Create CIVIC AID6 GKS representation."""
    clin_sig_prop = deepcopy(gks_therapeutic_proposition)
    clin_sig_prop["predicate"] = "hasClinicalSignificanceFor"
    clin_sig_prop["type"] = "VariantClinicalSignificanceProposition"
    clin_sig_prop.pop("objectTherapeutic")
    clin_sig_prop["objectCondition"] = clin_sig_prop.pop("conditionQualifier")

    params = {
        "id": "civic.aid:6",
        "description": "L858R is among the most common sensitizing EGFR mutations in NSCLC, and is assessed via DNA mutational analysis, including Sanger sequencing and next generation sequencing methods. Tyrosine kinase inhibitor afatinib is FDA approved as a first line systemic therapy in NSCLC with sensitizing EGFR mutation (civic.EID:2997).",
        "type": "Statement",
        "specifiedBy": gks_method,
        "proposition": clin_sig_prop,
        "direction": "supports",
        "strength": {
            "primaryCoding": {
                "system": "AMP/ASCO/CAP Guidelines, 2017",
                "code": "strong",
            },
        },
        "classification": {
            "name": "Tier I",
            "primaryCoding": {
                "system": "AMP/ASCO/CAP Guidelines, 2017",
                "code": "tier i",
            },
        },
        "hasEvidenceLines": [
            {
                "type": "EvidenceLine",
                "hasEvidenceItems": [gks_eid2997],
                "directionOfEvidenceProvided": "supports",
                "targetProposition": gks_therapeutic_proposition,
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "A",
                        "system": "AMP/ASCO/CAP Guidelines, 2017",
                    },
                },
            }
        ],
        "reportedIn": [
            "https://civicdb.org/links/assertion/6",
            gks_source592,
            {
                "type": "Document",
                "id": "civic.sid:592",
                "name": "Sequist et al., 2013",
                "title": "Phase III study of afatinib or cisplatin plus pemetrexed in patients with metastatic lung adenocarcinoma with EGFR mutations.",
                "pmid": "23816960",
                "urls": [
                    "https://civicdb.org/links/evidence/879",
                    "https://civicdb.org/links/source/592",
                    "http://www.ncbi.nlm.nih.gov/pubmed/23816960",
                ],
            },
            {
                "type": "Document",
                "id": "civic.sid:679",
                "name": "Wu et al., 2014",
                "title": "Afatinib versus cisplatin plus gemcitabine for first-line treatment of Asian patients with advanced non-small-cell lung cancer harbouring EGFR mutations (LUX-Lung 6): an open-label, randomised phase 3 trial.",
                "pmid": "24439929",
                "urls": [
                    "https://civicdb.org/links/evidence/982",
                    "https://civicdb.org/links/source/679",
                    "http://www.ncbi.nlm.nih.gov/pubmed/24439929",
                ],
            },
            {
                "type": "Document",
                "id": "civic.sid:594",
                "name": "Yang et al., 2012",
                "title": "Afatinib for patients with lung adenocarcinoma and epidermal growth factor receptor mutations (LUX-Lung 2): a phase 2 trial.",
                "pmid": "22452895",
                "urls": [
                    "https://civicdb.org/links/evidence/883",
                    "https://civicdb.org/links/source/594",
                    "http://www.ncbi.nlm.nih.gov/pubmed/22452895",
                ],
            },
            {
                "type": "Document",
                "id": "civic.sid:669",
                "name": "Hirano et al., 2015",
                "title": "In vitro modeling to determine mutation specificity of EGFR tyrosine kinase inhibitors against clinically relevant EGFR mutants in non-small-cell lung cancer.",
                "pmid": "26515464",
                "urls": [
                    "https://civicdb.org/links/evidence/968",
                    "https://civicdb.org/links/source/669",
                    "http://www.ncbi.nlm.nih.gov/pubmed/26515464",
                    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770737",
                ],
            },
            {
                "type": "Document",
                "id": "civic.sid:1525",
                "name": "Li et al., 2008",
                "title": "BIBW2992, an irreversible EGFR/HER2 inhibitor highly effective in preclinical lung cancer models.",
                "pmid": "18408761",
                "urls": [
                    "https://civicdb.org/links/evidence/2629",
                    "https://civicdb.org/links/source/1525",
                    "http://www.ncbi.nlm.nih.gov/pubmed/18408761",
                    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2748240",
                ],
            },
        ],
    }
    return VariantClinicalSignificanceStatement(**params)


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


def _ccv_method(method_type: str) -> dict:
    """Get CCV Method"""
    return {
        "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
        "reportedIn": {
            "id": "pmid:35101336",
            "name": "Horak et al., 2022, Genet Med.",
            "title": "Standards for the classification of pathogenicity of somatic variants in cancer (oncogenicity): Joint recommendations of Clinical Genome Resource (ClinGen), Cancer Genomics Consortium (CGC), and Variant Interpretation for Cancer Consortium (VICC)",
            "doi": "10.1016/j.gim.2022.01.001",
            "pmid": "35101336",
            "urls": [
                "https://doi.org/10.1016/j.gim.2022.01.001",
                "https://pubmed.ncbi.nlm.nih.gov/35101336/",
            ],
            "type": "Document",
        },
        "methodType": method_type,
        "type": "Method",
    }


@pytest.fixture(scope="module")
def gks_gid42():
    """Create test fixture for CIViC GID42 GKS representation."""
    return {
        "id": "civic.gid:42",
        "conceptType": "Gene",
        "name": "RET",
        "mappings": [
            {
                "coding": {
                    "id": "ncbigene:5979",
                    "code": "5979",
                    "system": "https://www.ncbi.nlm.nih.gov/gene/",
                },
                "relation": "exactMatch",
            },
        ],
        "extensions": [
            {
                "name": "description",
                "value": "RET mutations and the RET fusion RET-PTC lead to activation of this tyrosine kinase receptor and are associated with thyroid cancers. RET point mutations are the most common mutations identified in medullary thyroid cancer (MTC) with germline and somatic mutations in RET associated with hereditary and sporadic forms, respectively. The most common somatic form mutation is M918T (exon 16) and a variety of other mutations effecting exons 10, 11 and 15 have been described. The prognostic significance of these mutations have been hotly debated in the field, however, data suggests that some RET mutation may confer drug resistance. Highly selective and well-tolerated RET inhibitors, selpercatinib (LOXO-292) and pralsetinib (BLU-667), have been FDA approved recently for the treatment of RET fusion-positive non-small-cell lung cancer, RET fusion-positive thyroid cancer and RET-mutant medullary thyroid cancer.",
            },
            {
                "name": "aliases",
                "value": [
                    "CDHF12",
                    "CDHR16",
                    "HSCR1",
                    "MEN2A",
                    "MEN2B",
                    "MTC1",
                    "PTC",
                    "RET",
                    "RET-ELE1",
                ],
            },
        ],
    }


@pytest.fixture(scope="module")
def ret_m918t_vrs():
    return {
        "id": "ga4gh:VA.hEybNB_CeKflfFhT5AKOU5i1lgZPP-aS",
        "type": "Allele",
        "name": "RET M918T",
        "digest": "hEybNB_CeKflfFhT5AKOU5i1lgZPP-aS",
        "location": {
            "id": "ga4gh:SL.oIeqSfOEuqO7KNOPt8YUIa9vo1f6yMao",
            "type": "SequenceLocation",
            "digest": "oIeqSfOEuqO7KNOPt8YUIa9vo1f6yMao",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.jMu9-ItXSycQsm4hyABeW_UfSNRXRVnl",
            },
            "start": 917,
            "end": 918,
            "sequence": "M",
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "T"},
    }


@pytest.fixture(scope="module")
def civic_mpid113_cdna_vrs():
    return {
        "id": "ga4gh:VA.TZBjEPHhLRYxssQopcOQLWEBQrwzhH3T",
        "type": "Allele",
        "digest": "TZBjEPHhLRYxssQopcOQLWEBQrwzhH3T",
        "location": {
            "id": "ga4gh:SL.LD_QnJ8V1MR3stLat01acwyO4fWrUGco",
            "type": "SequenceLocation",
            "digest": "LD_QnJ8V1MR3stLat01acwyO4fWrUGco",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.jHlgYyFWJThVNL_o5UXEBwcQVNEPc62c",
            },
            "start": 2942,
            "end": 2943,
            "sequence": "T",
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "C"},
    }


@pytest.fixture(scope="module")
def civic_mpid113_genomic_vrs():
    return {
        "id": "ga4gh:VA.ON-Q17mJBYx3unmQ8GiqllzEphxR-Fie",
        "type": "Allele",
        "digest": "ON-Q17mJBYx3unmQ8GiqllzEphxR-Fie",
        "location": {
            "id": "ga4gh:SL.wIzpygPWdaZBkoKcIg461KaERW7XfyZS",
            "type": "SequenceLocation",
            "digest": "wIzpygPWdaZBkoKcIg461KaERW7XfyZS",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB",
            },
            "start": 43121967,
            "end": 43121968,
            "sequence": "T",
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "C"},
    }


@pytest.fixture(scope="module")
def civic_mpid113(ret_m918t_vrs, civic_mpid113_cdna_vrs, civic_mpid113_genomic_vrs):
    ret_m918t_vrs_copy = deepcopy(ret_m918t_vrs)
    hgvs_p = "NP_065681.1:p.Met918Thr"
    ret_m918t_vrs_copy["expressions"] = [
        {"syntax": "hgvs.p", "value": hgvs_p},
        {"syntax": "hgvs.p", "value": "ENSP00000347942.3:p.Met918Thr"},
    ]

    civic_mpid113_cdna_vrs_copy = deepcopy(civic_mpid113_cdna_vrs)
    hgvs_c = "NM_020975.4:c.2753T>C"
    civic_mpid113_cdna_vrs_copy["name"] = hgvs_c
    civic_mpid113_cdna_vrs_copy["expressions"] = [{"syntax": "hgvs.c", "value": hgvs_c}]

    civic_mpid113_genomic_vrs_copy = deepcopy(civic_mpid113_genomic_vrs)
    hgvs_g = "NC_000010.11:g.43121968T>C"
    civic_mpid113_genomic_vrs_copy["name"] = hgvs_g
    civic_mpid113_genomic_vrs_copy["expressions"] = [
        {"syntax": "hgvs.g", "value": hgvs_g},
        {"syntax": "hgvs.g", "value": "NC_000010.10:g.43617416T>C"},
    ]

    return {
        "id": "civic.mpid:113",
        "type": "CategoricalVariant",
        "description": "RET M918T is the most common somatically acquired mutation in medullary thyroid cancer (MTC). While there currently are no RET-specific inhibiting agents, promiscuous kinase inhibitors have seen some success in treating RET overactivity. Data suggests however, that the M918T mutation may lead to drug resistance, especially against the VEGFR-inhibitor motesanib. It has also been suggested that RET M918T leads to more aggressive MTC with a poorer prognosis.",
        "name": "RET M918T",
        "aliases": ["MET918THR"],
        "constraints": [
            {
                "type": "DefiningAlleleConstraint",
                "allele": ret_m918t_vrs_copy,
                "relations": [
                    {
                        "primaryCoding": {
                            "system": "ga4gh-gks-term:allele-relation",
                            "code": "liftover_to",
                        }
                    },
                    {
                        "primaryCoding": {
                            "system": "http://www.sequenceontology.org",
                            "code": "translation_of",
                        }
                    },
                ],
            }
        ],
        "mappings": [
            {
                "coding": {
                    "id": "dbsnp:rs74799832",
                    "code": "rs74799832",
                    "system": "https://www.ncbi.nlm.nih.gov/snp/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "clingen.allele:CA009082",
                    "code": "CA009082",
                    "system": "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "clinvar:13919",
                    "code": "13919",
                    "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "id": "civic.mpid:113",
                    "code": "113",
                    "system": "https://civicdb.org/links/molecular_profile/",
                },
                "relation": "exactMatch",
            },
            {
                "coding": {
                    "code": "113",
                    "id": "civic.vid:113",
                    "name": "M918T",
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
        "members": [
            civic_mpid113_cdna_vrs_copy,
            civic_mpid113_genomic_vrs_copy,
            ret_m918t_vrs_copy,
        ],
        "extensions": [
            {
                "name": "CIViC representative coordinate",
                "value": {
                    "chromosome": "10",
                    "start": 43617416,
                    "stop": 43617416,
                    "reference_bases": "T",
                    "variant_bases": "C",
                    "representative_transcript": "ENST00000355710.3",
                    "ensembl_version": 75,
                    "reference_build": "GRCh37",
                    "type": "coordinates",
                },
            },
            {
                "name": "CIViC Molecular Profile Score",
                "value": 139.0,
            },
            {
                "name": "hgvsDescriptions",
                "value": [
                    "ENST00000355710.3:c.2753T>C",
                    "NM_020975.4:c.2753T>C",
                    "NC_000010.10:g.43617416T>C",
                    "NP_065681.1:p.Met918Thr",
                    "NC_000010.11:g.43121968T>C",
                    "ENSP00000347942.3:p.Met918Thr",
                ],
            },
            {
                "name": "maneSelectTranscript",
                "value": "ENST00000355710.8:c.2753T>C",
            },
            {"name": "categoricalVariationType", "value": "ProteinSequenceConsequence"},
        ],
    }


@pytest.fixture(scope="module")
def civic_mpid113_undefined(civic_mpid113):
    """Create test fixture for undefined categorical variant"""

    mp = deepcopy(civic_mpid113)
    del mp["constraints"]
    del mp["members"]

    for ext in mp["extensions"]:
        if ext["name"] == "categoricalVariationType":
            ext["value"] = "Undefined"

    return mp


@pytest.fixture(scope="module")
def gks_aid202_proposition(gks_gid42, civic_mpid113):
    """Create test fixture forCIVIC AID6 proposition"""
    return {
        "type": "VariantOncogenicityProposition",
        "geneContextQualifier": gks_gid42,
        "objectTumorType": {
            "id": "civic.did:15",
            "conceptType": "Disease",
            "name": "Medullary Thyroid Carcinoma",
            "mappings": [
                {
                    "coding": {
                        "code": "DOID:3973",
                        "system": "https://disease-ontology.org/?id=",
                    },
                    "relation": "exactMatch",
                }
            ],
        },
        "alleleOriginQualifier": {
            "name": "somatic",
            "mappings": [
                {
                    "coding": {
                        "code": "SOMATIC",
                        "system": "https://civicdb.org",
                        "iris": [
                            "https://civic.readthedocs.io/en/latest/model/evidence/origin.html"
                        ],
                    },
                    "relation": "exactMatch",
                }
            ],
        },
        "predicate": "isOncogenicFor",
        "subjectVariant": civic_mpid113,
    }


@pytest.fixture(scope="module")
def gks_aid202(gks_aid202_proposition):
    """Create CIVIC AID6 GKS representation."""
    params = {
        "id": "civic.aid:202",
        "type": "Statement",
        "description": "Published sequencing studies have shown that RET mutations are very common in medullary thryoid carcinoma (MTC) and M918T is the most common specific variant, especially in the MEN2B clinical subtype of familial disease (civic.EID:78) but also in sporadic cases(civic.EID:12800). M918T mutations may predict worse outcomes (civic.EID:74). Biochemical and functional characterization demonstrates that the M918T mutation leads to functional activation of RET relative to wild-type through multiple complementary mechanisms, including increased ATP affinity (>10-fold) and complex stability, reduced conformational rigidity, and the promotion of ligand-independent dimerization and autophosphorylation (civic.EID:12805). Exogenous expression has been shown to induce transformation of Ba/F3 cells (civic.EID:11723), and drive colony formation in NIH3T3 cells (civic.EID:12709, OS2). RET M918T occurs in the region of the tyrosine kinase domain which is associated with multiple endocrine neoplasia type 2 B (OM1). RET M918T is predicted to be deleterious (CHASMplus score 0.314 > VECS gene-specific cutoff of 0.22, OP1). Eleven instances of the variant occur in cancerhotspots.org (V2): 6 Thyroid, 4 Adrenal Gland, 1 Breast (OP3). The variant is absent in gnomAD database (v4.1.0, OP4). Together these criteria indicate that M918T is likely oncogenic, with a score of 9.",
        "proposition": gks_aid202_proposition,
        "strength": {
            "primaryCoding": {
                "code": "likely",
                "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
            }
        },
        "classification": {
            "primaryCoding": {
                "code": "likely oncogenic",
                "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
            }
        },
        "reportedIn": [
            "https://civicdb.org/links/assertion/202",
            {
                "id": "civic.sid:44",
                "type": "Document",
                "name": "Elisei et al., 2008",
                "title": "Prognostic significance of somatic RET oncogene mutations in sporadic medullary thyroid cancer: a 10-year follow-up study.",
                "urls": [
                    "https://civicdb.org/links/evidence/74",
                    "https://civicdb.org/links/source/44",
                    "http://www.ncbi.nlm.nih.gov/pubmed/18073307",
                    "https://civicdb.org/links/evidence/12800",
                ],
                "pmid": "18073307",
            },
            {
                "id": "civic.sid:92",
                "type": "Document",
                "name": "Egawa et al., 1998",
                "title": "Genotype-phenotype correlation of patients with multiple endocrine neoplasia type 2 in Japan.",
                "urls": [
                    "https://civicdb.org/links/evidence/78",
                    "https://civicdb.org/links/source/92",
                    "http://www.ncbi.nlm.nih.gov/pubmed/9839497",
                ],
                "pmid": "9839497",
            },
            {
                "id": "civic.sid:5458",
                "type": "Document",
                "name": "Romei et al., 2018",
                "title": "RET mutation heterogeneity in primary advanced medullary thyroid cancers and their metastases.",
                "urls": [
                    "https://civicdb.org/links/evidence/12711",
                    "https://civicdb.org/links/source/5458",
                    "http://www.ncbi.nlm.nih.gov/pubmed/29515777",
                    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5839408",
                ],
                "pmid": "29515777",
            },
            {
                "id": "civic.sid:5519",
                "type": "Document",
                "name": "Gujral et al., 2006",
                "title": "Molecular mechanisms of RET receptor-mediated oncogenesis in multiple endocrine neoplasia 2B.",
                "urls": [
                    "https://civicdb.org/links/evidence/12805",
                    "https://civicdb.org/links/source/5519",
                    "http://www.ncbi.nlm.nih.gov/pubmed/17108110",
                ],
                "pmid": "17108110",
            },
            {
                "id": "civic.sid:4870",
                "type": "Document",
                "name": "Zhao et al., 2020",
                "title": "Identifying novel oncogenic RET mutations and characterising their sensitivity to RET-specific inhibitors.",
                "urls": [
                    "https://civicdb.org/links/evidence/11723",
                    "https://civicdb.org/links/source/4870",
                    "http://www.ncbi.nlm.nih.gov/pubmed/32284345",
                ],
                "pmid": "32284345",
            },
            {
                "id": "civic.sid:4953",
                "type": "Document",
                "name": "Ceccherini et al., 1997",
                "title": "Somatic in frame deletions not involving juxtamembranous cysteine residues strongly activate the RET proto-oncogene.",
                "urls": [
                    "https://civicdb.org/links/evidence/12709",
                    "https://civicdb.org/links/source/4953",
                    "http://www.ncbi.nlm.nih.gov/pubmed/9191060",
                ],
                "pmid": "9191060",
            },
        ],
        "direction": "supports",
        "specifiedBy": _ccv_method("guideline"),
        "hasEvidenceLines": [
            {
                "type": "EvidenceLine",
                "directionOfEvidenceProvided": "supports",
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "moderate",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "evidenceOutcome": {
                    "primaryCoding": {
                        "code": "OM1",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "scoreOfEvidenceProvided": 2,
                "specifiedBy": _ccv_method("functional_domain_location"),
            },
            {
                "type": "EvidenceLine",
                "directionOfEvidenceProvided": "supports",
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "strong",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "evidenceOutcome": {
                    "primaryCoding": {
                        "code": "OS2",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "scoreOfEvidenceProvided": 4,
                "specifiedBy": _ccv_method("functional_assay"),
            },
            {
                "type": "EvidenceLine",
                "directionOfEvidenceProvided": "supports",
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "supporting",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "evidenceOutcome": {
                    "primaryCoding": {
                        "code": "OP4",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "scoreOfEvidenceProvided": 1,
                "specifiedBy": _ccv_method("population_frequency"),
            },
            {
                "type": "EvidenceLine",
                "directionOfEvidenceProvided": "supports",
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "supporting",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "evidenceOutcome": {
                    "primaryCoding": {
                        "code": "OP1",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "scoreOfEvidenceProvided": 1,
                "specifiedBy": _ccv_method("computational_prediction"),
            },
            {
                "type": "EvidenceLine",
                "directionOfEvidenceProvided": "supports",
                "strengthOfEvidenceProvided": {
                    "primaryCoding": {
                        "code": "supporting",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "evidenceOutcome": {
                    "primaryCoding": {
                        "code": "OP3",
                        "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    }
                },
                "scoreOfEvidenceProvided": 1,
                "specifiedBy": _ccv_method("somatic_hotspot_recurrence"),
            },
        ],
    }
    return VariantOncogenicityStatement(**params)


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

    def test_valid(
        self,
        ret_m918t_vrs,
        civic_mpid113,
        civic_mpid113_cdna_vrs,
        civic_mpid113_genomic_vrs,
        monkeypatch,
    ):

        def normalize_side_effect(expr):
            if expr == "RET M918T":
                return Allele.model_validate(ret_m918t_vrs)

            if expr == "NM_020975.4:c.2753T>C":
                return Allele.model_validate(civic_mpid113_cdna_vrs)

            if expr in ["NC_000010.10:g.43617416T>C", "NC_000010.11:g.43121968T>C"]:
                return Allele.model_validate(civic_mpid113_genomic_vrs)

            if expr in [
                "ENST00000355710.3:c.2753T>C",
                "ENST00000355710.8:c.2753T>C",
            ]:
                return None

            raise AssertionError(f"Unexpected normalize query: {expr}")

        mp = civic.get_molecular_profile_by_id(113)

        variation_normalizer = Mock(spec=VariationNormalizerDataProxy)
        variation_normalizer.normalize.side_effect = normalize_side_effect
        variation_normalizer.normalize_molecular_profile.side_effect = (
            lambda molecular_profile: variation_normalizer.normalize("RET M918T")
        )

        monkeypatch.setattr(
            CivicGksMolecularProfile, "_variation_normalizer", variation_normalizer
        )
        gks_mp = CivicGksMolecularProfile(molecular_profile=mp)

        assert variation_normalizer.normalize.call_count == 6

        diff = DeepDiff(
            gks_mp.model_dump(exclude_none=True),
            civic_mpid113,
            ignore_order=True,
        )

        assert diff == {}

    def test_extensions_preserve_source_hgvs_descriptions(
        self, v600e_mp, mocked_normalizer
    ):
        """Source HGVS descriptions are retained even when not VRS expressions."""
        variant = v600e_mp.variants[0]

        with (
            patch.object(
                variant,
                "hgvs_expressions",
                new=["N/A", "XR_001744858.1:n.1823-3918T>A"],
            ),
            patch.object(
                variant,
                "mane_select_transcript",
                new=None,
            ),
        ):
            gks_mp = CivicGksMolecularProfile(v600e_mp)
            extensions = gks_mp.extensions
            assert extensions

            extension_values = {ext.name: ext.value for ext in extensions}
            assert extension_values["hgvsDescriptions"] == [
                "N/A",
                "XR_001744858.1:n.1823-3918T>A",
            ]
            assert "maneSelectTranscript" not in extension_values

    def test_expressions_deduplicate_mane_select(self, v600e_mp):
        """A MANE Select expression already in HGVS descriptions is emitted once."""
        variant = v600e_mp.variants[0]
        mane_select = "NM_005228.5:c.1799T>A"

        with (
            patch.object(variant, "hgvs_expressions", new=[mane_select]),
            patch.object(variant, "mane_select_transcript", new=mane_select),
        ):
            expressions = CivicGksMolecularProfile._get_expressions(variant)

        assert len(expressions) == 1
        assert expressions[0].value == mane_select
        assert expressions[0].extensions == [Extension(name="isManeSelect", value=True)]

    def test_na_clinvar_mapping(self, v600e_mp, mocked_normalizer):
        """Test that get_aliases_and_mappings method works as expected when no clinvar entry found"""
        variant = v600e_mp.variants[0]

        with patch.object(variant, "clinvar_entries", new=["N/A"]):
            gks_mp = CivicGksMolecularProfile(v600e_mp)
            mappings = gks_mp.mappings
            assert mappings
            assert not any(
                m.coding.system == "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
                for m in mappings
            )

    def test_build_constraints_gene_mutation(self, mocked_normalizer):
        mp = civic.get_molecular_profile_by_id(395)
        gks_mp = CivicGksMolecularProfile(mp)
        constraints = gks_mp.constraints
        assert constraints
        assert len(constraints) == 1
        assert gks_mp.members is None
        assert constraints[0].model_dump(exclude_none=True) == {
            "type": "FeatureContextConstraint",
            "featureContext": {
                "primaryCoding": {
                    "code": "673",
                    "id": "ncbigene:673",
                    "system": "https://www.ncbi.nlm.nih.gov/gene/",
                }
            },
        }

        categorical_variation_type = next(
            (
                ext.value
                for ext in gks_mp.extensions
                if ext.name == "categoricalVariationType"
            ),
            None,
        )
        assert categorical_variation_type == "FeatureContext"

    @patch.object(
        CivicGksGene,
        "get_mappings",
    )
    def test_build_constraints_no_gene_mappings(
        self,
        test_get_mappings,
        mocked_normalizer,
    ):
        test_get_mappings.return_value = []

        mp = civic.get_molecular_profile_by_id(395)

        with pytest.raises(
            CivicGksRecordError,
            match="Unable to retrieve mappings for gene 5",
        ):
            CivicGksMolecularProfile(mp)

    def test_build_constraints_normalization_failure(
        self, civic_mpid113_undefined, monkeypatch
    ):
        mp = civic.get_molecular_profile_by_id(113)

        variation_normalizer = Mock()
        variation_normalizer.normalize_molecular_profile.return_value = None

        monkeypatch.setattr(
            CivicGksMolecularProfile, "_variation_normalizer", variation_normalizer
        )
        gks_mp = CivicGksMolecularProfile(mp)
        diff = DeepDiff(
            gks_mp.model_dump(exclude_none=True),
            civic_mpid113_undefined,
            ignore_order=True,
        )
        assert diff == {}

    def test_build_constraints_allele(self, braf_v600e_vrs, mocked_normalizer):
        mp = civic.get_molecular_profile_by_id(12)
        normalized_allele = braf_v600e_vrs.model_dump(exclude_none=True)

        gks_mp = CivicGksMolecularProfile(mp)
        assert braf_v600e_vrs.model_dump(exclude_none=True) == normalized_allele

        constraints = gks_mp.constraints
        assert constraints
        assert len(constraints) == 1
        constraint = constraints[0].model_dump(exclude_none=True)
        protein_expressions = constraint["allele"].pop("expressions")
        assert constraint == {
            "type": "DefiningAlleleConstraint",
            "allele": normalized_allele,
            "relations": [
                {
                    "primaryCoding": {
                        "code": "liftover_to",
                        "system": "ga4gh-gks-term:allele-relation",
                    }
                },
                {
                    "primaryCoding": {
                        "code": "translation_of",
                        "system": "http://www.sequenceontology.org",
                    }
                },
            ],
        }
        assert protein_expressions
        assert {expression["syntax"] for expression in protein_expressions} == {
            "hgvs.p"
        }

    def test_build_members_does_not_mutate_normalized_variation(self, braf_v600e_vrs):
        """Member metadata is added to a copy of the normalizer result."""
        normalized_allele = braf_v600e_vrs.model_dump(exclude_none=True)
        expression = Expression(syntax=Syntax.HGVS_C, value="NM_004333.6:c.1799T>A")
        variation_normalizer = Mock()
        variation_normalizer.normalize.return_value = braf_v600e_vrs
        original_normalizer = CivicGksMolecularProfile._variation_normalizer

        try:
            CivicGksMolecularProfile.configure_variation_normalizer(
                variation_normalizer
            )
            members = CivicGksMolecularProfile._build_members(
                [expression], [Syntax.HGVS_C], None
            )
        finally:
            CivicGksMolecularProfile._variation_normalizer = original_normalizer

        assert braf_v600e_vrs.model_dump(exclude_none=True) == normalized_allele
        assert len(members) == 1
        assert members[0].root.name == expression.value
        assert members[0].root.expressions == [expression]

    def test_build_constraints_copy_number_change(
        self,
        braf_amplification_vrs,
        monkeypatch,
    ):
        mp = civic.get_molecular_profile_by_id(1243)

        variation_normalizer = Mock()
        variation_normalizer.normalize_molecular_profile.return_value = (
            braf_amplification_vrs
        )

        monkeypatch.setattr(
            CivicGksMolecularProfile, "_variation_normalizer", variation_normalizer
        )
        gks_mp = CivicGksMolecularProfile(mp)
        constraints = gks_mp.constraints
        assert constraints
        assert len(constraints) == 2
        constraints_dict = [c.model_dump(exclude_none=True) for c in constraints]
        diff = DeepDiff(
            constraints_dict,
            [
                {
                    "type": "CopyChangeConstraint",
                    "copyChange": braf_amplification_vrs.copyChange,
                },
                {
                    "type": "DefiningLocationConstraint",
                    "location": braf_amplification_vrs.location.model_dump(
                        exclude_none=True
                    ),
                    "relations": [
                        {
                            "primaryCoding": {
                                "code": "liftover_to",
                                "system": "ga4gh-gks-term:allele-relation",
                            }
                        }
                    ],
                    "matchCharacteristic": {
                        "primaryCoding": {
                            "code": "is_within",
                            "system": "ga4gh-gks-term:location-match",
                        }
                    },
                },
            ],
            ignore_order=True,
        )
        assert diff == {}

        categorical_variation_type = next(
            (
                ext.value
                for ext in gks_mp.extensions
                if ext.name == "categoricalVariationType"
            ),
            None,
        )
        assert categorical_variation_type == "CategoricalCnv"

    def test_build_constraints_unsupported_vrs_type(
        self,
        monkeypatch,
    ):
        # Need to pick non- Gene Mutation MP ID
        # Mocked value is dummy value purely for test purposes
        mp = civic.get_molecular_profile_by_id(1243)
        variation_normalizer = Mock()
        variation_normalizer.normalize_molecular_profile.return_value = CopyNumberCount(
            copies=1, location=iriReference("#/location/1")
        )

        monkeypatch.setattr(
            CivicGksMolecularProfile, "_variation_normalizer", variation_normalizer
        )
        with pytest.raises(
            CivicGksRecordError,
            match="Unsupported VRS variation type returned by Variation Normalizer. mpid=1243, type='CopyNumberCount'",
        ):
            CivicGksMolecularProfile(mp)

    def test_no_representative_coordinates(self):
        """Test that empty representative coordinates do not get an extension"""
        gks_mp = CivicGksMolecularProfile(civic.get_molecular_profile_by_id(2261))
        assert gks_mp.extensions
        rep_coord_ext = next(
            (
                ext
                for ext in gks_mp.extensions
                if ext.name == "CIViC representative coordinate"
            ),
            None,
        )
        assert rep_coord_ext is None


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


class TestCivicGksClinSigAssertion(object):
    """Test that CivicGksClinSigAssertion works as expected"""

    @patch("civicpy.exports.civic_gks_record.CivicGksMolecularProfile")
    def test_valid_single_therapy(self, test_mp, aid6, gks_aid6, gks_mpid33):
        """Test that single therapy works as expected"""
        test_mp.return_value = CategoricalVariant.model_validate(gks_mpid33)
        record = CivicGksClinSigAssertion(aid6)
        assert isinstance(record, VariantClinicalSignificanceStatement)
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

    def test_valid_combination_therapy(self, aid7, mocked_normalizer):
        """Test that combination therapy works as expected"""
        record = CivicGksClinSigAssertion(aid7)
        assert isinstance(record, VariantClinicalSignificanceStatement)
        assert len(record.hasEvidenceLines) == 1
        assert len(record.hasEvidenceLines[0].hasEvidenceItems) == 4
        therapy = record.hasEvidenceLines[0].targetProposition.objectTherapeutic.root
        assert isinstance(therapy, TherapyGroup)
        assert therapy.membershipOperator == "AND"
        assert len(therapy.therapies) == 2
        therapy_ids = {t.id for t in therapy.therapies}
        assert therapy_ids == {"civic.tid:19", "civic.tid:22"}

    @patch.object(civic.Assertion, "is_valid_for_gks_json")
    @patch.object(civic.Assertion, "evidence_items")
    @patch.object(civic.FusionVariant, "hgvs_expressions", create=True)
    @patch.object(civic.FusionVariant, "mane_select_transcript", create=True)
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
    @patch("civicpy.exports.civic_gks_record.CivicGksMolecularProfile")
    def test_valid_substitution_therapy(
        self,
        test_mp,
        test_coordinates,
        test_clinvar_entries,
        test_allele_registry_id,
        test_mane_select_transcript,
        test_hgvs_expressions,
        test_evidence_items,
        test_is_valid_for_gks_json,
        aid19,
        civic_mpid113,
    ):
        """Test that substitution therapy works as expected"""
        test_mp.return_value = CategoricalVariant.model_validate(civic_mpid113)
        test_coordinates.return_value = None
        test_clinvar_entries.return_value = []
        test_allele_registry_id.return_value = None
        test_is_valid_for_gks_json.return_value = True
        test_evidence_items.return_value = []
        test_hgvs_expressions.return_value = None
        test_mane_select_transcript.return_value = None
        record = CivicGksClinSigAssertion(aid19)
        assert isinstance(record, VariantClinicalSignificanceStatement)
        assert len(record.hasEvidenceLines) == 1
        therapy = record.hasEvidenceLines[0].targetProposition.objectTherapeutic.root
        assert isinstance(therapy, TherapyGroup)
        assert therapy.membershipOperator == "OR"
        assert len(therapy.therapies) == 2
        therapy_ids = {t.id for t in therapy.therapies}
        assert therapy_ids == {"civic.tid:5", "civic.tid:20"}

    def test_valid_prognostic(self, aid20, mocked_normalizer):
        """Test that valid prognostic assertion works as expected"""
        record = CivicGksClinSigAssertion(aid20)
        assert isinstance(record, VariantClinicalSignificanceStatement)
        assert len(record.hasEvidenceLines) == 1
        assert len(record.hasEvidenceLines[0].hasEvidenceItems) == 6
        assert (
            record.hasEvidenceLines[0].targetProposition.predicate
            == "associatedWithWorseOutcomeFor"
        )
        assert record.strength.primaryCoding.code.root == "strong"
        assert record.classification.primaryCoding.code.root == "tier i"

    @patch.object(civic.Assertion, "evidence_items", new_callable=PropertyMock)
    @patch.object(civic.Evidence, "is_valid_for_gks_json")
    def test_citations(
        self, test_is_valid_for_gks_json, test_evidence_items, aid20, mocked_normalizer
    ):
        """Test that citations extension is working correctly for EIDs that are not valid for GKS"""
        test_evidence_items.return_value = [civic.get_evidence_by_id(11881)]
        test_is_valid_for_gks_json.return_value = False

        record = CivicGksClinSigAssertion(aid20)
        assert len(record.hasEvidenceLines) == 1
        assert record.hasEvidenceLines[0].hasEvidenceItems is None

        reported_in = []
        for r in record.reportedIn:
            if isinstance(r, iriReference):
                reported_in.append(r.root)
            else:
                reported_in.append(r.model_dump(exclude_none=True))

        assert reported_in == [
            "https://civicdb.org/links/assertion/20",
            {
                "type": "Document",
                "id": "civic.sid:4914",
                "name": "Grimwade et al., 1998",
                "title": "The importance of diagnostic cytogenetics on outcome in AML: analysis of 1,612 patients entered into the MRC AML 10 trial. The Medical Research Council Adult and Children's Leukaemia Working Parties.",
                "pmid": "9746770",
                "urls": [
                    "https://civicdb.org/links/evidence/11881",
                    "https://civicdb.org/links/source/4914",
                    "http://www.ncbi.nlm.nih.gov/pubmed/9746770",
                ],
            },
        ]


class TestCivicGksDiagnosticAssertion(object):
    """Test that CivicGksDiagnosticAssertion works as expected"""

    def test_valid(
        self,
        mocked_normalizer,
        aid9,
        aid93,
        gks_aid93_object_condition,
        aid115,
        gks_aid115_object_condition,
    ):
        """Test that valid diagnostic assertion works as expected"""
        record = CivicGksClinSigAssertion(aid9)
        assert isinstance(record, VariantClinicalSignificanceStatement)
        assert len(record.hasEvidenceLines) == 1
        assert len(record.hasEvidenceLines[0].hasEvidenceItems) == 2
        assert (
            record.hasEvidenceLines[0].targetProposition.predicate
            == "isDiagnosticInclusionCriterionFor"
        )
        assert record.strength.primaryCoding.code.root == "potential"
        assert record.classification.primaryCoding.code.root == "tier ii"
        assert (
            record.hasEvidenceLines[
                0
            ].strengthOfEvidenceProvided.primaryCoding.code.root
            == "C"
        )

        # Single phenotype (complex condition set)
        record = CivicGksClinSigAssertion(aid93)
        assert isinstance(record, VariantClinicalSignificanceStatement)
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
        record = CivicGksClinSigAssertion(aid115)
        assert isinstance(record, VariantClinicalSignificanceStatement)
        record_object_condition = record.proposition.objectCondition
        assert isinstance(record_object_condition, Condition)
        assert isinstance(record_object_condition.root, ConditionSet)
        diff = DeepDiff(
            record_object_condition.model_dump(exclude_none=True),
            gks_aid115_object_condition,
            ignore_order=True,
        )
        assert diff == {}

    def test_invalid(self, aid117):
        """Test that unsupported assertion types raise exceptions"""

        with pytest.raises(
            CivicGksRecordError,
            match=re.escape(
                "Assertion type must be one of ['PREDICTIVE', 'PROGNOSTIC', 'DIAGNOSTIC']"
            ),
        ):
            CivicGksClinSigAssertion(aid117)


class TestCivicGksOncogenicAssertion(object):
    """Test that CivicGksOncogenicAssertion works as expected"""

    @patch("civicpy.exports.civic_gks_record.CivicGksMolecularProfile")
    def test_valid(
        self,
        test_mp,
        civic_mpid113,
        aid202,
        gks_aid202,
    ):
        """Test that valid oncogenic assertions works as expected"""

        def evidence_key(item: dict) -> str:
            return item["evidenceOutcome"]["primaryCoding"]["code"]

        test_mp.return_value = CategoricalVariant.model_validate(civic_mpid113)
        record = CivicGksOncogenicAssertion(aid202, approval=None)

        assert isinstance(record, VariantOncogenicityStatement)

        actual = record.model_dump(exclude_none=True)
        expected = gks_aid202.model_dump(exclude_none=True)

        assert set(actual.keys()) == set(expected.keys())

        # Split out due to large record
        for key in expected:
            if key == "hasEvidenceLines":
                actual_evidence = actual[key]
                expected_evidence = expected[key]

                assert len(actual_evidence) == len(expected_evidence), (
                    f"Mismatch in hasEvidenceLines length: "
                    f"actual={len(actual_evidence)}, expected={len(expected_evidence)}"
                )

                actual_by_code = {evidence_key(item): item for item in actual_evidence}
                expected_by_code = {
                    evidence_key(item): item for item in expected_evidence
                }

                assert set(actual_by_code) == set(expected_by_code), (
                    "Mismatch in hasEvidence evidenceOutcome.primaryCoding.code values"
                )

                for code in expected_by_code:
                    diff = DeepDiff(
                        actual_by_code[code],
                        expected_by_code[code],
                        ignore_order=True,
                    )

                    assert diff == {}, (
                        "Mismatch in hasEvidence item with "
                        f"evidenceOutcome.primaryCoding.code={code}"
                    )

                continue

            diff = DeepDiff(
                actual[key],
                expected[key],
                ignore_order=True,
            )

            assert diff == {}, f"Mismatch in key: {key}"

    def test_invalid(self, aid6):
        """Test that unsupported assertion types raise exceptions"""

        with pytest.raises(
            CivicGksRecordError,
            match=re.escape("Assertion type must be one of ['ONCOGENIC']"),
        ):
            CivicGksOncogenicAssertion(aid6)


class TestCivicGksRecord(object):
    """Test that GKS Record helper functions work correctly"""

    @patch("civicpy.exports.civic_gks_record.VariationNormalizerRESTDataProxy")
    def test_configured_variation_normalizer_is_shared(self, mock_rest_normalizer):
        original_normalizer = CivicGksMolecularProfile._variation_normalizer
        variation_normalizer = Mock(spec=VariationNormalizerDataProxy)

        try:
            CivicGksMolecularProfile.configure_variation_normalizer(
                variation_normalizer
            )

            assert (
                CivicGksMolecularProfile._get_variation_normalizer()
                is variation_normalizer
            )
            assert (
                CivicGksMolecularProfile._get_variation_normalizer()
                is variation_normalizer
            )
            assert "_variation_normalizer" not in CivicGksMolecularProfile.model_fields
            mock_rest_normalizer.assert_not_called()
        finally:
            CivicGksMolecularProfile._variation_normalizer = original_normalizer

    @patch("civicpy.exports.civic_gks_record.VariationNormalizerRESTDataProxy")
    def test_default_variation_normalizer_is_initialized_lazily(
        self, mock_rest_normalizer
    ):
        original_normalizer = CivicGksMolecularProfile._variation_normalizer
        variation_normalizer = Mock(spec=VariationNormalizerDataProxy)
        mock_rest_normalizer.return_value = variation_normalizer

        try:
            CivicGksMolecularProfile._variation_normalizer = None

            assert (
                CivicGksMolecularProfile._get_variation_normalizer()
                is variation_normalizer
            )
            assert (
                CivicGksMolecularProfile._get_variation_normalizer()
                is variation_normalizer
            )
            mock_rest_normalizer.assert_called_once_with()
        finally:
            CivicGksMolecularProfile._variation_normalizer = original_normalizer

    def test_unsupported_assertion_type(self, mocked_normalizer):
        """Test that unsupported assertion types raise NotImplementedError"""

        with pytest.raises(
            NotImplementedError,
            match=r"Assertion type PREDISPOSING is not currently supported",
        ):
            create_gks_record_from_assertion(civic.get_assertion_by_id(17))

    @patch("civicpy.exports.civic_gks_record.CivicGksClinSigAssertion")
    def test_factory_preserves_positional_approval(self, mock_gks_assertion, aid6):
        """The factory retains its existing positional approval argument."""
        approval = Mock()

        create_gks_record_from_assertion(aid6, approval)

        mock_gks_assertion.assert_called_once_with(
            aid6,
            approval=approval,
        )

    @pytest.mark.parametrize(
        (
            "civic_assertion_fixture_name",
            "submission_type_filter",
            "should_raise_error",
        ),
        (
            [
                "aid202",
                ClinVarSubmissionType.ONCOGENICITY,
                False,
            ],
            [
                "aid202",
                ClinVarSubmissionType.CLINICAL_IMPACT,
                True,
            ],
            [
                "aid9",
                ClinVarSubmissionType.CLINICAL_IMPACT,
                False,
            ],
            [
                "aid9",
                ClinVarSubmissionType.ONCOGENICITY,
                True,
            ],
            [
                "aid20",
                ClinVarSubmissionType.CLINICAL_IMPACT,
                False,
            ],
            [
                "aid20",
                ClinVarSubmissionType.ONCOGENICITY,
                True,
            ],
            [
                "aid6",
                ClinVarSubmissionType.CLINICAL_IMPACT,
                False,
            ],
            [
                "aid6",
                ClinVarSubmissionType.ONCOGENICITY,
                True,
            ],
        ),
    )
    def test_create_gks_record_from_assertion_filter(
        self,
        request,
        civic_assertion_fixture_name,
        submission_type_filter,
        should_raise_error,
        mocked_normalizer,
    ):
        """Test that create_gks_record_from_assertion works correctly when submission filter is applied"""
        civic_aid = request.getfixturevalue(civic_assertion_fixture_name)
        if should_raise_error:
            with pytest.raises(
                NotImplementedError,
                match=rf"Assertion type {civic_aid.assertion_type} is not supported for ClinVar submission type {submission_type_filter.value}",
            ):
                create_gks_record_from_assertion(
                    civic_aid,
                    submission_type_filter=submission_type_filter,
                )
        else:
            assert create_gks_record_from_assertion(
                civic_aid,
                submission_type_filter=submission_type_filter,
            )

    def test_clinvar_accession_ext(self, mocked_normalizer):
        a = civic.get_assertion_by_id(193)
        record = create_gks_record_from_assertion(
            a,
            approval=a.approvals[0],
        )
        assert isinstance(record, VariantClinicalSignificanceStatement)
        assert [ext.model_dump(exclude_none=True) for ext in record.extensions] == [
            {"name": "clinvar_accession", "value": "SCV007542591"}
        ]

    def test_assertion_invalid(self, aid117, mocked_normalizer):
        """Test that invalid assertion raise exceptions"""

        with pytest.raises(
            CivicGksRecordError,
            match=re.escape(
                "Assertion type must be one of ['PREDICTIVE', 'PROGNOSTIC', 'DIAGNOSTIC']"
            ),
        ):
            CivicGksClinSigAssertion(aid117)
