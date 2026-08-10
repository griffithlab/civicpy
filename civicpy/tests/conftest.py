from copy import deepcopy
from typing import Any
from unittest.mock import Mock

import pytest
from ga4gh.vrs.models import Allele, CopyNumberChange
from ga4gh.va_spec.aac_2017 import VariantClinicalSignificanceStatement
from ga4gh.va_spec.base import Statement
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement

from civicpy.civic import Assertion, get_assertion_by_id, get_molecular_profile_by_id
from civicpy.exports.civic_gks_record import CivicGksMolecularProfile
from civicpy.exports.variation_normalizer import VariationNormalizerDataProxy


@pytest.fixture(scope="module")
def braf_v600e_vrs():
    params = {
        "name": "BRAF V600E",
        "id": "ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
        "type": "Allele",
        "digest": "j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
        "location": {
            "id": "ga4gh:SL.t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
            "type": "SequenceLocation",
            "digest": "t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
            },
            "start": 599,
            "end": 600,
            "sequence": "V",
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
    }
    return Allele(**params)


@pytest.fixture(scope="session")
def braf_amplification_vrs():
    params = {
        "name": "BRAF Amplification",
        "id": "ga4gh:CX.h7unj-f_djER28-h2Q6Prvo3C90O4d3M",
        "type": "CopyNumberChange",
        "digest": "h7unj-f_djER28-h2Q6Prvo3C90O4d3M",
        "location": {
            "id": "ga4gh:SL.0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
            "type": "SequenceLocation",
            "digest": "0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            },
            "start": 140713327,
            "end": 140924929,
        },
        "copyChange": "high-level gain",
    }
    return CopyNumberChange(**params)


@pytest.fixture(scope="module")
def v600e_mp():
    return get_molecular_profile_by_id(12)


@pytest.fixture(scope="module")
def aid6() -> Assertion:
    """Return the predictive Assertion shared by GKS export tests."""
    return get_assertion_by_id(6)


@pytest.fixture(scope="module")
def aid202() -> Assertion:
    """Return the oncogenic Assertion shared by GKS export tests."""
    return get_assertion_by_id(202)


@pytest.fixture(scope="module")
def mocked_normalizer(braf_v600e_vrs):
    """Provide a normalizer mock for GKS tests without a live service."""
    variation_normalizer = Mock(spec=VariationNormalizerDataProxy)
    variation_normalizer.normalize.return_value = None
    variation_normalizer.normalize_molecular_profile.return_value = braf_v600e_vrs
    previous_normalizer = CivicGksMolecularProfile._variation_normalizer
    CivicGksMolecularProfile.configure_variation_normalizer(variation_normalizer)
    yield variation_normalizer
    CivicGksMolecularProfile._variation_normalizer = previous_normalizer


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
                            "name": "variantTypes",
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
                "name": "representativeVariantCoordinates",
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
                "name": "molecularProfileScore",
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
                            "name": "variantTypes",
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
                "name": "representativeVariantCoordinates",
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
                "name": "molecularProfileScore",
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


def _find_gks_objects(value: Any, object_type: str) -> list[dict[str, Any]]:
    """Find typed GKS objects nested in a fixture representation."""
    matches: list[dict[str, Any]] = []
    if isinstance(value, dict):
        if value.get("type") == object_type and "id" in value:
            matches.append(deepcopy(value))
        for nested_value in value.values():
            matches.extend(_find_gks_objects(nested_value, object_type))
    elif isinstance(value, list):
        for nested_value in value:
            matches.extend(_find_gks_objects(nested_value, object_type))
    return matches


@pytest.fixture(scope="module")
def gks_bundle_sequence_references(braf_v600e_vrs) -> dict[str, Any]:
    """Return the sequence reference from the shared normalized allele."""
    allele = braf_v600e_vrs.model_dump(exclude_none=True)
    sequence_reference = deepcopy(allele["location"]["sequenceReference"])
    return {sequence_reference["refgetAccession"]: sequence_reference}


@pytest.fixture(scope="module")
def gks_bundle_locations(braf_v600e_vrs) -> dict[str, Any]:
    """Reference the sequence used by the shared normalized allele location."""
    allele = braf_v600e_vrs.model_dump(exclude_none=True)
    location = deepcopy(allele["location"])
    refget_accession = location["sequenceReference"]["refgetAccession"]
    location["sequenceReference"] = f"#/sequenceReference/{refget_accession}"
    return {location["id"]: location}


@pytest.fixture(scope="module")
def gks_bundle_variants(
    braf_v600e_vrs,
    gks_mpid33,
    civic_mpid113,
) -> dict[str, Any]:
    """Return the VRS object under each associated CIViC variant ID."""
    allele = braf_v600e_vrs.model_dump(exclude_none=True)
    allele["location"] = f"#/location/{allele['location']['id']}"
    allele_33 = deepcopy(allele)
    allele_33["expressions"] = deepcopy(
        gks_mpid33["constraints"][0]["allele"]["expressions"]
    )
    allele_33["expressions"].append(
        {"syntax": "hgvs.p", "value": "ENSP00000275493.2:p.Leu858Arg"}
    )
    allele["expressions"] = deepcopy(
        civic_mpid113["constraints"][0]["allele"]["expressions"]
    )
    return {
        "civic.vid:33": {"protein": {allele_33["id"]: allele_33}},
        "civic.vid:113": {"protein": {allele["id"]: allele}},
    }


@pytest.fixture(scope="module")
def gks_bundle_features(gks_gid19, gks_gid42) -> dict[str, Any]:
    """Return genes already defined by the shared GKS record fixtures."""
    return {gene["id"]: deepcopy(gene) for gene in (gks_gid19, gks_gid42)}


@pytest.fixture(scope="module")
def gks_bundle_molecular_profiles() -> dict[str, Any]:
    """Return the expected ``molecularProfile`` bundle collection."""
    return {
        "civic.mpid:113": {
            "id": "civic.mpid:113",
            "type": "CategoricalVariant",
            "name": "RET M918T",
            "description": "RET M918T is the most common somatically acquired "
            "mutation in medullary thyroid cancer (MTC). While "
            "there currently are no RET-specific inhibiting "
            "agents, promiscuous kinase inhibitors have seen "
            "some success in treating RET overactivity. Data "
            "suggests however, that the M918T mutation may lead "
            "to drug resistance, especially against the "
            "VEGFR-inhibitor motesanib. It has also been "
            "suggested that RET M918T leads to more aggressive "
            "MTC with a poorer prognosis.",
            "aliases": ["MET918THR"],
            "extensions": [
                {"name": "molecularProfileScore", "value": 139.0},
                {
                    "name": "hgvsDescriptions",
                    "value": [
                        "NM_020975.4:c.2753T>C",
                        "NP_065681.1:p.Met918Thr",
                        "ENST00000355710.3:c.2753T>C",
                        "NC_000010.10:g.43617416T>C",
                        "NC_000010.11:g.43121968T>C",
                        "ENSP00000347942.3:p.Met918Thr",
                    ],
                },
                {
                    "name": "maneSelectTranscript",
                    "value": "ENST00000355710.8:c.2753T>C",
                },
                {
                    "name": "representativeVariantCoordinates",
                    "value": {
                        "chromosome": "10",
                        "start": 43617416,
                        "stop": 43617416,
                        "reference_bases": "T",
                        "variant_bases": "C",
                        "ensembl_version": 75,
                        "representative_transcript": "ENST00000355710.3",
                        "reference_build": "GRCh37",
                        "type": "coordinates",
                    },
                },
                {
                    "name": "categoricalVariationType",
                    "value": "ProteinSequenceConsequence",
                },
            ],
            "members": [
                "#/variant/civic.vid:113/protein/ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L"
            ],
            "constraints": [
                {
                    "type": "DefiningAlleleConstraint",
                    "allele": "#/variant/civic.vid:113/protein/ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
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
                    "coding": "#/molecularProfile/civic.mpid:113",
                    "relation": "exactMatch",
                },
                {
                    "coding": {
                        "id": "civic.vid:113",
                        "extensions": [
                            {"name": "subtype", "value": "gene_variant"},
                            {
                                "name": "variantTypes",
                                "value": [
                                    {
                                        "coding": {
                                            "id": "civic.variant_type:47",
                                            "name": "Missense Variant",
                                            "system": "http://www.sequenceontology.org/browser/current_svn/term/",
                                            "code": "SO:0001583",
                                        },
                                        "relation": "exactMatch",
                                    }
                                ],
                            },
                        ],
                        "name": "M918T",
                        "system": "https://civicdb.org/links/variant/",
                        "code": "113",
                    },
                    "relation": "exactMatch",
                },
                {
                    "coding": {
                        "id": "clingen.allele:CA009082",
                        "system": "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=",
                        "code": "CA009082",
                    },
                    "relation": "relatedMatch",
                },
                {
                    "coding": {
                        "id": "clinvar:13919",
                        "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                        "code": "13919",
                    },
                    "relation": "relatedMatch",
                },
                {
                    "coding": {
                        "id": "dbsnp:rs74799832",
                        "system": "https://www.ncbi.nlm.nih.gov/snp/",
                        "code": "rs74799832",
                    },
                    "relation": "relatedMatch",
                },
            ],
        },
        "civic.mpid:33": {
            "id": "civic.mpid:33",
            "type": "CategoricalVariant",
            "name": "EGFR L858R",
            "description": "EGFR L858R has long been recognized as a "
            "functionally significant mutation in cancer, and is "
            "one of the most prevalent single mutations in lung "
            "cancer. Best described in non-small cell lung "
            "cancer (NSCLC), the mutation seems to confer "
            "sensitivity to first and second generation TKI's "
            "like gefitinib and neratinib. NSCLC patients with "
            "this mutation treated with TKI's show increased "
            "overall and progression-free survival, as compared "
            "to chemotherapy alone. Third generation TKI's are "
            "currently in clinical trials that specifically "
            "focus on mutant forms of EGFR, a few of which have "
            "shown efficacy in treating patients that failed to "
            "respond to earlier generation TKI therapies.",
            "aliases": ["LEU858ARG", "L813R", "LEU813ARG"],
            "extensions": [
                {"name": "molecularProfileScore", "value": 379.0},
                {
                    "name": "hgvsDescriptions",
                    "value": [
                        "NC_000007.13:g.55259515T>G",
                        "NM_005228.4:c.2573T>G",
                        "ENST00000275493.2:c.2573T>G",
                        "NP_005219.2:p.Leu858Arg",
                        "NC_000007.14:g.55191822T>G",
                        "ENSP00000275493.2:p.Leu858Arg",
                    ],
                },
                {
                    "name": "maneSelectTranscript",
                    "value": "ENST00000275493.7:c.2573T>G",
                },
                {
                    "name": "representativeVariantCoordinates",
                    "value": {
                        "chromosome": "7",
                        "start": 55259515,
                        "stop": 55259515,
                        "reference_bases": "T",
                        "variant_bases": "G",
                        "ensembl_version": 75,
                        "representative_transcript": "ENST00000275493.2",
                        "reference_build": "GRCh37",
                        "type": "coordinates",
                    },
                },
                {
                    "name": "categoricalVariationType",
                    "value": "ProteinSequenceConsequence",
                },
            ],
            "members": [
                "#/variant/civic.vid:33/protein/ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L"
            ],
            "constraints": [
                {
                    "type": "DefiningAlleleConstraint",
                    "allele": "#/variant/civic.vid:33/protein/ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
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
                    "coding": "#/molecularProfile/civic.mpid:33",
                    "relation": "exactMatch",
                },
                {
                    "coding": {
                        "id": "civic.vid:33",
                        "extensions": [
                            {"name": "subtype", "value": "gene_variant"},
                            {
                                "name": "variantTypes",
                                "value": [
                                    {
                                        "coding": {
                                            "id": "civic.variant_type:47",
                                            "name": "Missense Variant",
                                            "system": "http://www.sequenceontology.org/browser/current_svn/term/",
                                            "code": "SO:0001583",
                                        },
                                        "relation": "exactMatch",
                                    }
                                ],
                            },
                        ],
                        "name": "L858R",
                        "system": "https://civicdb.org/links/variant/",
                        "code": "33",
                    },
                    "relation": "exactMatch",
                },
                {
                    "coding": {
                        "id": "clingen.allele:CA126713",
                        "system": "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=",
                        "code": "CA126713",
                    },
                    "relation": "relatedMatch",
                },
                {
                    "coding": {
                        "id": "clinvar:16609",
                        "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                        "code": "16609",
                    },
                    "relation": "relatedMatch",
                },
                {
                    "coding": {
                        "id": "clinvar:376282",
                        "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                        "code": "376282",
                    },
                    "relation": "relatedMatch",
                },
                {
                    "coding": {
                        "id": "clinvar:376280",
                        "system": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                        "code": "376280",
                    },
                    "relation": "relatedMatch",
                },
                {
                    "coding": {
                        "id": "dbsnp:rs121434568",
                        "system": "https://www.ncbi.nlm.nih.gov/snp/",
                        "code": "rs121434568",
                    },
                    "relation": "relatedMatch",
                },
            ],
        },
    }


@pytest.fixture(scope="module")
def gks_bundle_diseases(gks_did8) -> dict[str, Any]:
    """Return diseases shared by the record and bundle expectations."""
    diseases = {
        "civic.did:15": {
            "id": "civic.did:15",
            "conceptType": "Disease",
            "name": "Medullary Thyroid Carcinoma",
            "mappings": [
                {
                    "coding": {
                        "system": "https://disease-ontology.org/?id=",
                        "code": "DOID:3973",
                    },
                    "relation": "exactMatch",
                }
            ],
        },
        "civic.did:30": {
            "id": "civic.did:30",
            "conceptType": "Disease",
            "name": "Lung Adenocarcinoma",
            "mappings": [
                {
                    "coding": {
                        "system": "https://disease-ontology.org/?id=",
                        "code": "DOID:3910",
                    },
                    "relation": "exactMatch",
                }
            ],
        },
    }
    diseases[gks_did8["id"]] = deepcopy(gks_did8)
    return diseases


@pytest.fixture(scope="module")
def gks_bundle_condition_sets() -> dict[str, Any]:
    """Return the expected ``conditionSet`` bundle collection."""
    return {}


@pytest.fixture(scope="module")
def gks_bundle_therapies(gks_tid146) -> dict[str, Any]:
    """Return therapies already defined by the shared GKS record fixtures."""
    return {gks_tid146["id"]: deepcopy(gks_tid146)}


@pytest.fixture(scope="module")
def gks_bundle_therapy_groups() -> dict[str, Any]:
    """Return the expected ``therapyGroup`` bundle collection."""
    return {}


@pytest.fixture(scope="module")
def gks_bundle_variant_origins() -> dict[str, Any]:
    """Return the expected ``variantOrigin`` bundle collection."""
    return {
        "civic.variantOrigin:SOMATIC": {
            "id": "civic.variantOrigin:SOMATIC",
            "name": "somatic",
            "mappings": [
                {
                    "coding": {
                        "system": "https://civicdb.org",
                        "code": "SOMATIC",
                        "iris": [
                            "https://civic.readthedocs.io/en/latest/model/evidence/origin.html"
                        ],
                    },
                    "relation": "exactMatch",
                }
            ],
        }
    }


@pytest.fixture(scope="module")
def gks_bundle_sources(gks_method, gks_aid6, gks_aid202) -> dict[str, Any]:
    """Collect documents from the shared dereferenced GKS fixtures."""
    documents: dict[str, Any] = {}
    for value in (
        gks_method,
        gks_aid6.model_dump(exclude_none=True),
        gks_aid202.model_dump(exclude_none=True),
    ):
        for document in _find_gks_objects(value, "Document"):
            documents[document["id"]] = document
    return documents


@pytest.fixture(scope="module")
def gks_bundle_methods(gks_method) -> dict[str, Any]:
    """Reference the document nested in the shared method fixture."""
    method = deepcopy(gks_method)
    method["reportedIn"] = f"#/source/{method['reportedIn']['id']}"
    return {method["id"]: method}


@pytest.fixture(scope="module")
def gks_bundle_organizations() -> dict[str, Any]:
    """Return the expected ``organization`` bundle collection."""
    return {}


@pytest.fixture(scope="module")
def gks_bundle_propositions() -> dict[str, Any]:
    """Return the expected ``proposition`` bundle collection."""
    return {
        "civic.proposition:-AKWXtNluL_XZYk5cDaaV7bKw6fKlPmD": {
            "id": "civic.proposition:-AKWXtNluL_XZYk5cDaaV7bKw6fKlPmD",
            "type": "VariantClinicalSignificanceProposition",
            "subjectVariant": "#/molecularProfile/civic.mpid:33",
            "geneContextQualifier": "#/feature/civic.gid:19",
            "alleleOriginQualifier": "#/variantOrigin/civic.variantOrigin:SOMATIC",
            "predicate": "hasClinicalSignificanceFor",
            "objectCondition": "#/disease/civic.did:8",
        },
        "civic.proposition:lGNyTBSVq9ncomifdlwtOERYq7ZM37FX": {
            "id": "civic.proposition:lGNyTBSVq9ncomifdlwtOERYq7ZM37FX",
            "type": "VariantOncogenicityProposition",
            "subjectVariant": "#/molecularProfile/civic.mpid:113",
            "geneContextQualifier": "#/feature/civic.gid:42",
            "alleleOriginQualifier": "#/variantOrigin/civic.variantOrigin:SOMATIC",
            "predicate": "isOncogenicFor",
            "objectTumorType": "#/disease/civic.did:15",
        },
        "civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu": {
            "id": "civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu",
            "type": "VariantTherapeuticResponseProposition",
            "subjectVariant": "#/molecularProfile/civic.mpid:33",
            "geneContextQualifier": "#/feature/civic.gid:19",
            "alleleOriginQualifier": "#/variantOrigin/civic.variantOrigin:SOMATIC",
            "predicate": "predictsSensitivityTo",
            "objectTherapeutic": "#/therapy/civic.tid:146",
            "conditionQualifier": "#/disease/civic.did:8",
        },
        "civic.proposition:nqvFeEaF3J52FxWjzgqxOoafka3s50pY": {
            "id": "civic.proposition:nqvFeEaF3J52FxWjzgqxOoafka3s50pY",
            "type": "VariantTherapeuticResponseProposition",
            "subjectVariant": "#/molecularProfile/civic.mpid:33",
            "geneContextQualifier": "#/feature/civic.gid:19",
            "alleleOriginQualifier": "#/variantOrigin/civic.variantOrigin:SOMATIC",
            "predicate": "predictsSensitivityTo",
            "objectTherapeutic": "#/therapy/civic.tid:146",
            "conditionQualifier": "#/disease/civic.did:30",
        },
    }


@pytest.fixture(scope="module")
def gks_bundle_statement_objects() -> dict[str, Any]:
    """Return the expected ``statement`` bundle collection."""
    return {
        "civic.aid:202": {
            "id": "civic.aid:202",
            "type": "Statement",
            "description": "Published sequencing studies have shown that RET "
            "mutations are very common in medullary thryoid "
            "carcinoma (MTC) and M918T is the most common "
            "specific variant, especially in the MEN2B clinical "
            "subtype of familial disease (civic.EID:78) but also "
            "in sporadic cases(civic.EID:12800). M918T mutations "
            "may predict worse outcomes (civic.EID:74). "
            "Biochemical and functional characterization "
            "demonstrates that the M918T mutation leads to "
            "functional activation of RET relative to wild-type "
            "through multiple complementary mechanisms, "
            "including increased ATP affinity (>10-fold) and "
            "complex stability, reduced conformational rigidity, "
            "and the promotion of ligand-independent "
            "dimerization and autophosphorylation "
            "(civic.EID:12805). Exogenous expression has been "
            "shown to induce transformation of Ba/F3 cells "
            "(civic.EID:11723), and drive colony formation in "
            "NIH3T3 cells (civic.EID:12709, OS2). RET M918T "
            "occurs in the region of the tyrosine kinase domain "
            "which is associated with multiple endocrine "
            "neoplasia type 2 B (OM1). RET M918T is predicted to "
            "be deleterious (CHASMplus score 0.314 > VECS "
            "gene-specific cutoff of 0.22, OP1). Eleven "
            "instances of the variant occur in "
            "cancerhotspots.org (V2): 6 Thyroid, 4 Adrenal "
            "Gland, 1 Breast (OP3). The variant is absent in "
            "gnomAD database (v4.1.0, OP4). Together these "
            "criteria indicate that M918T is likely oncogenic, "
            "with a score of 9.",
            "specifiedBy": {
                "type": "Method",
                "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                "methodType": "guideline",
                "reportedIn": "#/source/pmid:35101336",
            },
            "reportedIn": [
                "https://civicdb.org/links/assertion/202",
                "#/source/civic.sid:44",
                "#/source/civic.sid:92",
                "#/source/civic.sid:5458",
                "#/source/civic.sid:5519",
                "#/source/civic.sid:4870",
                "#/source/civic.sid:4953",
            ],
            "proposition": "#/proposition/civic.proposition:lGNyTBSVq9ncomifdlwtOERYq7ZM37FX",
            "direction": "supports",
            "strength": {
                "primaryCoding": {
                    "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    "code": "likely",
                }
            },
            "classification": {
                "primaryCoding": {
                    "system": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                    "code": "likely oncogenic",
                }
            },
            "hasEvidenceLines": [
                {
                    "type": "EvidenceLine",
                    "specifiedBy": {
                        "type": "Method",
                        "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                        "methodType": "functional_domain_location",
                        "reportedIn": "#/source/pmid:35101336",
                    },
                    "directionOfEvidenceProvided": "supports",
                    "strengthOfEvidenceProvided": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "moderate",
                        }
                    },
                    "scoreOfEvidenceProvided": 2,
                    "evidenceOutcome": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "OM1",
                        }
                    },
                },
                {
                    "type": "EvidenceLine",
                    "specifiedBy": {
                        "type": "Method",
                        "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                        "methodType": "functional_assay",
                        "reportedIn": "#/source/pmid:35101336",
                    },
                    "directionOfEvidenceProvided": "supports",
                    "strengthOfEvidenceProvided": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "strong",
                        }
                    },
                    "scoreOfEvidenceProvided": 4,
                    "evidenceOutcome": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "OS2",
                        }
                    },
                },
                {
                    "type": "EvidenceLine",
                    "specifiedBy": {
                        "type": "Method",
                        "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                        "methodType": "population_frequency",
                        "reportedIn": "#/source/pmid:35101336",
                    },
                    "directionOfEvidenceProvided": "supports",
                    "strengthOfEvidenceProvided": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "supporting",
                        }
                    },
                    "scoreOfEvidenceProvided": 1,
                    "evidenceOutcome": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "OP4",
                        }
                    },
                },
                {
                    "type": "EvidenceLine",
                    "specifiedBy": {
                        "type": "Method",
                        "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                        "methodType": "computational_prediction",
                        "reportedIn": "#/source/pmid:35101336",
                    },
                    "directionOfEvidenceProvided": "supports",
                    "strengthOfEvidenceProvided": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "supporting",
                        }
                    },
                    "scoreOfEvidenceProvided": 1,
                    "evidenceOutcome": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "OP1",
                        }
                    },
                },
                {
                    "type": "EvidenceLine",
                    "specifiedBy": {
                        "type": "Method",
                        "name": "ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022",
                        "methodType": "somatic_hotspot_recurrence",
                        "reportedIn": "#/source/pmid:35101336",
                    },
                    "directionOfEvidenceProvided": "supports",
                    "strengthOfEvidenceProvided": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "supporting",
                        }
                    },
                    "scoreOfEvidenceProvided": 1,
                    "evidenceOutcome": {
                        "primaryCoding": {
                            "system": "ClinGen/CGC/VICC "
                            "Guidelines "
                            "for "
                            "Oncogenicity, "
                            "2022",
                            "code": "OP3",
                        }
                    },
                },
            ],
        },
        "civic.aid:6": {
            "id": "civic.aid:6",
            "type": "Statement",
            "description": "L858R is among the most common sensitizing EGFR "
            "mutations in NSCLC, and is assessed via DNA "
            "mutational analysis, including Sanger sequencing and "
            "next generation sequencing methods. Tyrosine kinase "
            "inhibitor afatinib is FDA approved as a first line "
            "systemic therapy in NSCLC with sensitizing EGFR "
            "mutation (civic.EID:2997).",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": [
                "https://civicdb.org/links/assertion/6",
                "#/source/civic.sid:1725",
                "#/source/civic.sid:592",
                "#/source/civic.sid:679",
                "#/source/civic.sid:594",
                "#/source/civic.sid:669",
                "#/source/civic.sid:1525",
            ],
            "proposition": "#/proposition/civic.proposition:-AKWXtNluL_XZYk5cDaaV7bKw6fKlPmD",
            "direction": "supports",
            "strength": {
                "primaryCoding": {
                    "system": "AMP/ASCO/CAP Guidelines, 2017",
                    "code": "strong",
                }
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
                    "targetProposition": "#/proposition/civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu",
                    "hasEvidenceItems": [
                        "#/evidence/civic.eid:2997",
                        "#/evidence/civic.eid:879",
                        "#/evidence/civic.eid:982",
                        "#/evidence/civic.eid:883",
                        "#/evidence/civic.eid:968",
                        "#/evidence/civic.eid:2629",
                    ],
                    "directionOfEvidenceProvided": "supports",
                    "strengthOfEvidenceProvided": {
                        "primaryCoding": {
                            "system": "AMP/ASCO/CAP Guidelines, 2017",
                            "code": "A",
                        }
                    },
                }
            ],
        },
        "civic.eid:2629": {
            "id": "civic.eid:2629",
            "type": "Statement",
            "description": "In an in vitro study using NCI-H1666 cells "
            "(wildtype EGFR) and NCI-H3255 cells (EGFR-L858R), "
            "inhibition of cell growth was used as an assay to "
            "determine sensitivity to irreversible tyrosine "
            "kinase inhibitor (TKI) drugs. Cells with an EGFR "
            "L858R mutation demonstrated an improved response "
            "to afatinib (IC50: 0.7nM vs. 60nM) compared to "
            "wildtype EGFR cells.",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": ["#/source/civic.sid:1525"],
            "proposition": "#/proposition/civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu",
            "direction": "supports",
            "strength": {
                "name": "Preclinical evidence",
                "primaryCoding": {
                    "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                    "code": "D",
                },
                "mappings": [
                    {
                        "coding": {
                            "name": "preclinical evidence",
                            "system": "https://go.osu.edu/evidence-codes",
                            "code": "e000009",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        },
        "civic.eid:2997": {
            "id": "civic.eid:2997",
            "type": "Statement",
            "description": "Afatinib, an irreversible inhibitor of the ErbB "
            "family of tyrosine kinases has been approved in "
            "the US for the first-line treatment of patients "
            "with metastatic non-small-cell lung cancer (NSCLC) "
            "who have tumours with EGFR exon 19 deletions or "
            "exon 21 (L858R) substitution mutations as detected "
            "by a US FDA-approved test",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": ["#/source/civic.sid:1725"],
            "proposition": "#/proposition/civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu",
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
                            "name": "authoritative evidence",
                            "system": "https://go.osu.edu/evidence-codes",
                            "code": "e000001",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        },
        "civic.eid:879": {
            "id": "civic.eid:879",
            "type": "Statement",
            "description": "A phase III clinical trial (NCT00949650) found that "
            "median progression free survival among patients "
            "with exon 19 deletions or L858R EGFR mutations (n = "
            "308) was 13.6 months for afatinib and 6.9 months "
            "for chemotherapy (HR, 0.47; 95% CI, 0.34 to 0.65; P "
            "= 0.001).",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": ["#/source/civic.sid:592"],
            "proposition": "#/proposition/civic.proposition:nqvFeEaF3J52FxWjzgqxOoafka3s50pY",
            "direction": "supports",
            "strength": {
                "name": "Clinical evidence",
                "primaryCoding": {
                    "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                    "code": "B",
                },
                "mappings": [
                    {
                        "coding": {
                            "name": "clinical cohort evidence",
                            "system": "https://go.osu.edu/evidence-codes",
                            "code": "e000005",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        },
        "civic.eid:883": {
            "id": "civic.eid:883",
            "type": "Statement",
            "description": "In a phase 2 study of patients with lung "
            "adenocarcinoma (stage IIIb with pleural effusion or "
            "stage IV) and EGFR mutations, treated with afatinib "
            "were assessed by objective response. 129 patients "
            "were treated with afatinib. 66% of the 106 patients "
            "with two common activating EGFR mutations (deletion "
            "19 or L858R) had an objective response compared to "
            "39% of 23 patients with less common mutations.",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": ["#/source/civic.sid:594"],
            "proposition": "#/proposition/civic.proposition:nqvFeEaF3J52FxWjzgqxOoafka3s50pY",
            "direction": "supports",
            "strength": {
                "name": "Clinical evidence",
                "primaryCoding": {
                    "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                    "code": "B",
                },
                "mappings": [
                    {
                        "coding": {
                            "name": "clinical cohort evidence",
                            "system": "https://go.osu.edu/evidence-codes",
                            "code": "e000005",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        },
        "civic.eid:968": {
            "id": "civic.eid:968",
            "type": "Statement",
            "description": "Cells harboring L858R were sensitive to afatinib. "
            "This study performed drug response assays using "
            "five human NSCLC cell lines with various "
            "combinations of EGFR mutations. In order to "
            "directly compare the sensitivity of multiple EGFR "
            "mutations to EGFR-TKIs the authors also generated "
            "multiple EGFR transduced Ba/F3 stable cell lines "
            "and evaluated sensitivity to EGFR-TKIs by MTS "
            "assay.",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": ["#/source/civic.sid:669"],
            "proposition": "#/proposition/civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu",
            "direction": "supports",
            "strength": {
                "name": "Preclinical evidence",
                "primaryCoding": {
                    "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                    "code": "D",
                },
                "mappings": [
                    {
                        "coding": {
                            "name": "preclinical evidence",
                            "system": "https://go.osu.edu/evidence-codes",
                            "code": "e000009",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        },
        "civic.eid:982": {
            "id": "civic.eid:982",
            "type": "Statement",
            "description": "Afatinib is an irreversible covalent inhibitor of "
            "EGFR (second generation). This Phase III clinical "
            "trial (LUX-Lung 6; NCT01121393) was performed in "
            "Asian patients with EGFR mutant advanced NSCLC. 364 "
            "eligible patients with EGFR mutations were assigned "
            "to afatinib (n=242) or gemcitabine and cisplatin "
            "(n=122) treatment. The trial observed significantly "
            "longer median progression-free survival with "
            "afatinib vs. gemcitabine and cisplatin treatment "
            "(11.0 vs. 5.6 months). Afatinib/Chemotherapy group "
            "compositions: 51.2/50.8 % del 19; 38/37.7 % "
            "Leu858Arg; 10.8/11.5 % Uncommon.",
            "specifiedBy": "#/method/civic.method:2019",
            "reportedIn": ["#/source/civic.sid:679"],
            "proposition": "#/proposition/civic.proposition:nqvFeEaF3J52FxWjzgqxOoafka3s50pY",
            "direction": "supports",
            "strength": {
                "name": "Clinical evidence",
                "primaryCoding": {
                    "system": "https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                    "code": "B",
                },
                "mappings": [
                    {
                        "coding": {
                            "name": "clinical cohort evidence",
                            "system": "https://go.osu.edu/evidence-codes",
                            "code": "e000005",
                        },
                        "relation": "exactMatch",
                    }
                ],
            },
        },
    }


@pytest.fixture(scope="module")
def gks_bundle_metadata() -> dict[str, Any]:
    """Return the expected ``metadata`` bundle collection."""
    return {
        "va_spec_python_version": "test",
        "created_at": "2026-08-03",
        "bundle_format": "civic-gks-bundle",
        "bundle_format_version": "0.1.0",
        "statistics": {
            "collections": {
                "sequenceReference": {"count": 1},
                "location": {"count": 1, "types": {"SequenceLocation": 1}},
                "variant": {"count": 2, "types": {"Allele": 2}},
                "feature": {"count": 2},
                "molecularProfile": {"count": 2},
                "disease": {"count": 3},
                "phenotype": {"count": 0},
                "conditionSet": {"count": 0},
                "therapy": {"count": 1},
                "therapyGroup": {"count": 0},
                "variantOrigin": {"count": 1},
                "source": {"count": 14},
                "method": {"count": 1},
                "organization": {"count": 0},
                "proposition": {
                    "count": 4,
                    "types": {
                        "VariantClinicalSignificanceProposition": 1,
                        "VariantOncogenicityProposition": 1,
                        "VariantTherapeuticResponseProposition": 2,
                    },
                },
                "evidence": {"count": 6},
                "assertion": {"count": 2},
            }
        },
    }


@pytest.fixture(scope="module")
def gks_bundle_expected(
    gks_bundle_sequence_references,
    gks_bundle_locations,
    gks_bundle_variants,
    gks_bundle_features,
    gks_bundle_molecular_profiles,
    gks_bundle_diseases,
    gks_bundle_condition_sets,
    gks_bundle_therapies,
    gks_bundle_therapy_groups,
    gks_bundle_variant_origins,
    gks_bundle_sources,
    gks_bundle_methods,
    gks_bundle_organizations,
    gks_bundle_propositions,
    gks_bundle_statement_objects,
    gks_bundle_metadata,
) -> dict[str, Any]:
    """Compose the complete expected bundle from its collection fixtures."""
    return {
        "sequenceReference": gks_bundle_sequence_references,
        "location": gks_bundle_locations,
        "variant": gks_bundle_variants,
        "feature": gks_bundle_features,
        "molecularProfile": gks_bundle_molecular_profiles,
        "disease": gks_bundle_diseases,
        "phenotype": {},
        "conditionSet": gks_bundle_condition_sets,
        "therapy": gks_bundle_therapies,
        "therapyGroup": gks_bundle_therapy_groups,
        "variantOrigin": gks_bundle_variant_origins,
        "source": gks_bundle_sources,
        "method": gks_bundle_methods,
        "organization": gks_bundle_organizations,
        "proposition": gks_bundle_propositions,
        "evidence": {
            identifier: statement
            for identifier, statement in gks_bundle_statement_objects.items()
            if identifier.startswith("civic.eid:")
        },
        "assertion": {
            identifier: statement
            for identifier, statement in gks_bundle_statement_objects.items()
            if identifier.startswith("civic.aid:")
        },
        "metadata": gks_bundle_metadata,
        "failed_assertion_ids": [],
        "errors": [],
    }
