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
from ga4gh.va_spec.base import Condition, ConditionSet, TherapyGroup
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
def gks_contributions():
    return [
        {
            "type": "Contribution",
            "contributor": {
                "id": "civic.organization:1",
                "type": "Agent",
                "name": "CIViC",
                "description": "The CIViC Organization (formerly “The McDonnell Genome Institute” CIViC organization) comprises the founders, developers, editors, curators, and administrators who build and maintain the knowledgebase, based at Washington University in St. Louis. This group is dedicated to ensuring that high-quality cancer variant interpretations are broadly accessible for precision oncology. One of their main roles is evaluating and synthesizing crowdsourced community contributions into formal clinical Assertions. Once these Assertions meet the strict criteria of the CIViC standard operating procedure, core approval members approve them for 1-star submission to the ClinVar CIViC organization.",
                "extensions": [{"name": "isApprovedVcep", "value": False}],
            },
            "activityType": "approval.last_reviewed",
            "date": "2026-04-24",
        }
    ]


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
                if ext.name == "representativeVariantCoordinates"
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

    @pytest.mark.parametrize("criterion", ["OS2_moderate", "OS2_Moderate"])
    def test_va_spec_accepts_criterion_strength_case(self, criterion: str) -> None:
        """Accept lowercase and title-case criterion strength suffixes."""
        assertion = Mock(
            assertion_direction="SUPPORTS",
            clingen_codes=[Mock(code=criterion)],
        )

        evidence_lines = CivicGksOncogenicAssertion.get_evidence_lines(assertion)

        assert [
            line.model_dump(exclude_none=True)["evidenceOutcome"]["primaryCoding"][
                "code"
            ]
            for line in evidence_lines
        ] == ["OS2_moderate"]

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
            {"name": "clinvarAccession", "value": "SCV007542591"}
        ]
        contribution_extensions = record.contributions[0].extensions
        assert contribution_extensions
        assert [
            extension.model_dump(exclude_none=True)
            for extension in contribution_extensions
        ] == [{"name": "clinvarAccession", "value": "SCV007542591"}]
        assert all(
            isinstance(extension.name, str) and isinstance(extension.value, str)
            for extension in (*record.extensions, *contribution_extensions)
        )

    def test_assertion_invalid(self, aid117, mocked_normalizer):
        """Test that invalid assertion raise exceptions"""

        with pytest.raises(
            CivicGksRecordError,
            match=re.escape(
                "Assertion type must be one of ['PREDICTIVE', 'PROGNOSTIC', 'DIAGNOSTIC']"
            ),
        ):
            CivicGksClinSigAssertion(aid117)
