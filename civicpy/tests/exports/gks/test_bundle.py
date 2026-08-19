"""Tests for dereferenced and referenced CIViC GKS output documents."""

from datetime import date
from unittest.mock import Mock

import pytest
from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.base import Statement
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
from ga4gh.vrs.models import Allele, SequenceReference
from pydantic import ValidationError

from civicpy import civic
from civicpy.exports.civic_gks_record import create_gks_record_from_assertion
from civicpy.exports.gks.models import (
    GksAssertionError,
    GksOutputMetadata,
)
from civicpy.exports.gks.bundle import GksBundle, build_gks_bundle


class TestCivicGksBundleOutput:
    """Test the optional, reference-linked GKS Bundle Format export."""

    def test_schema_describes_concrete_bundle_objects(self) -> None:
        """Expose upstream GKS models instead of arbitrary JSON objects."""
        schema = GksBundle.model_json_schema()
        properties = schema["properties"]

        collection_key_patterns = {
            "sequenceReference": r"^SQ\.[A-Za-z0-9_-]+$",
            "location": r"^ga4gh:SL\.[A-Za-z0-9_-]+$",
            "variant": r"^civic\.vid:[0-9]+$",
            "feature": r"^civic\.gid:[0-9]+$",
            "molecularProfile": r"^civic\.mpid:[0-9]+$",
            "disease": r"^civic\.did:[0-9]+$",
            "phenotype": r"^civic\.phenotype:[0-9]+$",
            "conditionSet": r"^civic\.conditionSet:[A-Za-z0-9_-]+$",
            "therapy": r"^civic\.tid:[0-9]+$",
            "therapyGroup": r"^civic\.therapyGroup:[A-Za-z0-9_-]+$",
            "variantOrigin": r"^civic\.variantOrigin:[A-Za-z0-9_-]+$",
            "source": r"^(civic\.sid|pmid):[0-9]+$",
            "method": r"^civic\.method:[A-Za-z0-9_-]+$",
            "organization": r"^civic\.organization:[A-Za-z0-9_-]+$",
            "proposition": r"^civic\.proposition:[A-Za-z0-9_-]+$",
            "evidence": r"^civic\.eid:[0-9]+$",
            "assertion": r"^civic\.aid:[0-9]+$",
        }
        for collection, key_pattern in collection_key_patterns.items():
            collection_schema = properties[collection]
            assert collection_schema["additionalProperties"] is False
            assert key_pattern in collection_schema["patternProperties"]

        assert (
            len(
                properties["assertion"]["patternProperties"][
                    collection_key_patterns["assertion"]
                ]["anyOf"]
            )
            == 2
        )
        variant_schema = properties["variant"]["patternProperties"][
            collection_key_patterns["variant"]
        ]
        assert variant_schema["additionalProperties"] is False
        assert set(variant_schema["properties"]) == {
            "protein",
            "coding",
            "genomic",
            "unclassified",
        }
        for representation_schema in variant_schema["properties"].values():
            assert len(representation_schema["additionalProperties"]["anyOf"]) == 2
        assert {
            "Adjacency",
            "CisPhasedBlock",
            "CopyNumberCount",
            "DerivativeMolecule",
            "StudyResult",
            "Terminus",
            "Variation",
        }.isdisjoint(schema["$defs"])

        definition_references: set[str] = set()

        def collect_definition_references(value: object) -> None:
            if isinstance(value, dict):
                reference = value.get("$ref")
                if isinstance(reference, str) and reference.startswith("#/$defs/"):
                    definition_references.add(reference.rsplit("/", 1)[-1])
                for nested_value in value.values():
                    collect_definition_references(nested_value)
            elif isinstance(value, list):
                for nested_value in value:
                    collect_definition_references(nested_value)

        collect_definition_references(schema)
        assert definition_references <= schema["$defs"].keys()

        public_schema = GksBundle.model_json_schema()
        public_definition_references: set[str] = set()
        definition_references = public_definition_references
        collect_definition_references(public_schema)
        assert public_definition_references <= public_schema["$defs"].keys()
        assert "Allele" not in public_schema["$defs"]
        assert "GksAllele" not in public_schema["$defs"]
        assert "VariantOncogenicityStatement" not in public_schema["$defs"]
        assert "GksVariantOncogenicityStatement" not in public_schema["$defs"]
        assert "MappableConcept" not in public_schema["$defs"]
        assert "LengthExpression" not in public_schema["$defs"]
        assert "CopyChangeConstraint" not in public_schema["$defs"]
        assert "DiagnosticPredicate" not in public_schema["$defs"]

        external_references: set[str] = set()

        def collect_external_references(value: object) -> None:
            if isinstance(value, dict):
                reference = value.get("$ref")
                if isinstance(reference, str) and reference.startswith(
                    "https://w3id.org/ga4gh/schema/"
                ):
                    external_references.add(reference)
                for nested_value in value.values():
                    collect_external_references(nested_value)
            elif isinstance(value, list):
                for nested_value in value:
                    collect_external_references(nested_value)

        collect_external_references(public_schema)
        assert Allele.schema_id() in external_references
        assert MappableConcept.schema_id() in external_references
        assert SequenceReference.schema_id() in external_references
        assert CategoricalVariant.schema_id() in external_references
        assert Statement.schema_id() in external_references
        assert VariantOncogenicityStatement.schema_id() in external_references

        sequence_reference_schema = public_schema["properties"]["sequenceReference"][
            "patternProperties"
        ][collection_key_patterns["sequenceReference"]]
        assert sequence_reference_schema == {"$ref": SequenceReference.schema_id()}

        feature_schema = public_schema["properties"]["feature"]["patternProperties"][
            collection_key_patterns["feature"]
        ]
        assert feature_schema == {
            "allOf": [
                {
                    "$ref": MappableConcept.schema_id()
                },
                {
                    "properties": {"conceptType": {"const": "Gene"}},
                    "required": ["conceptType"],
                },
            ]
        }

    def test_builds_empty_bundle_with_errors(self) -> None:
        """Retain errors and zero counts when no Statements are available."""
        error = GksAssertionError(
            assertion_id=1,
            message="Unsupported value: #/this/is/not/a/bundle/reference",
        )

        bundle = build_gks_bundle(
            [],
            GksOutputMetadata(va_spec_python_version="test", created_at="2026-08-03"),
            [error],
        )

        assert bundle.failed_assertion_ids == [1]
        assert bundle.errors == [error]
        assert all(
            statistics.count == 0
            for statistics in bundle.metadata.statistics.collections.values()
        )

    def test_builds_with_default_metadata_and_errors(self) -> None:
        """Provide convenient defaults for optional bundle context."""
        bundle = build_gks_bundle([])

        assert bundle.metadata.created_at == date.today().isoformat()
        assert bundle.errors == []
        assert bundle.failed_assertion_ids == []

    def test_builds_real_assertions_across_bundle_collections(
        self, mocked_normalizer: Mock
    ) -> None:
        """Exercise conversion and bundle extraction without a full JSON snapshot."""
        records = [
            create_gks_record_from_assertion(civic.get_assertion_by_id(assertion_id))
            for assertion_id in (6, 202)
        ]

        bundle = build_gks_bundle(records)

        assert set(bundle.assertion) == {"civic.aid:6", "civic.aid:202"}
        assert bundle.evidence
        assert bundle.molecularProfile
        assert bundle.variant
        assert bundle.feature
        assert bundle.disease
        assert bundle.source
        assert bundle.method
        assert bundle.proposition

        statements = [*bundle.assertion.values(), *bundle.evidence.values()]
        assert all(
            statement["proposition"].startswith("#/proposition/")
            for statement in statements
        )
        assert all(
            constraint["allele"].startswith("#/variant/")
            for profile in bundle.molecularProfile.values()
            for constraint in profile.get("constraints", [])
        )

        statistics = bundle.metadata.statistics.collections
        for collection_name in (
            "assertion",
            "evidence",
            "molecularProfile",
            "feature",
            "disease",
            "source",
            "method",
            "proposition",
        ):
            assert statistics[collection_name].count == len(
                getattr(bundle, collection_name)
            )

    def test_duplicate_object_retains_first_without_orphans(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        first_record = Mock()
        first_record.model_dump.return_value = {
            "id": "civic.aid:1",
            "type": "Statement",
            "subject": {
                "id": "civic.mpid:1",
                "type": "CategoricalVariant",
                "mappings": [{"coding": {"id": "civic.vid:1"}}],
                "constraints": [{"allele": {"id": "ga4gh:VA.first", "type": "Allele"}}],
            },
        }
        second_record = Mock()
        second_record.model_dump.return_value = {
            "id": "civic.aid:2",
            "type": "Statement",
            "subject": {
                "id": "civic.mpid:1",
                "type": "CategoricalVariant",
                "mappings": [{"coding": {"id": "civic.vid:1"}}],
                "constraints": [
                    {"allele": {"id": "ga4gh:VA.discarded", "type": "Allele"}}
                ],
            },
        }

        bundle = build_gks_bundle(
            [second_record, first_record],
            GksOutputMetadata(va_spec_python_version="test", created_at="2026-08-03"),
            [],
        )

        assert set(bundle.molecularProfile) == {"civic.mpid:1"}
        assert set(bundle.variant) == {"civic.vid:1"}
        assert "retaining the first bundle object" in caplog.text

    def test_keys_vrs_variant_by_civic_vid(self) -> None:
        """Store one defining VRS object per VID and preserve other members."""
        record = Mock()
        record.model_dump.return_value = {
            "id": "civic.aid:1",
            "type": "Statement",
            "subject": {
                "id": "civic.mpid:1",
                "type": "CategoricalVariant",
                "mappings": [{"coding": {"id": "civic.vid:42"}}],
                "constraints": [
                    {
                        "allele": {
                            "id": "ga4gh:VA.allele",
                            "type": "Allele",
                        }
                    }
                ],
                "members": [
                    {
                        "id": "ga4gh:VA.member",
                        "type": "Allele",
                        "expressions": [{"syntax": "hgvs.c", "value": "c.1T>C"}],
                    },
                    {
                        "id": "ga4gh:VA.allele",
                        "type": "Allele",
                        "expressions": [{"syntax": "hgvs.p", "value": "p.V1A"}],
                    },
                ],
            },
        }

        bundle = build_gks_bundle(
            [record],
            GksOutputMetadata(va_spec_python_version="test", created_at="2026-08-03"),
            [],
        )

        assert bundle.variant == {
            "civic.vid:42": {
                "protein": {
                    "ga4gh:VA.allele": {
                        "id": "ga4gh:VA.allele",
                        "type": "Allele",
                    }
                },
                "coding": {
                    "ga4gh:VA.member": {
                        "id": "ga4gh:VA.member",
                        "type": "Allele",
                        "expressions": [{"syntax": "hgvs.c", "value": "c.1T>C"}],
                    }
                },
            }
        }
        assert bundle.molecularProfile["civic.mpid:1"]["constraints"] == [
            {"allele": "#/variant/civic.vid:42/protein/ga4gh:VA.allele"}
        ]
        assert bundle.molecularProfile["civic.mpid:1"]["members"] == [
            "#/variant/civic.vid:42/coding/ga4gh:VA.member",
            "#/variant/civic.vid:42/protein/ga4gh:VA.allele",
        ]
        assert bundle.metadata.statistics.collections["variant"].types == {"Allele": 2}

    def test_references_condition_and_therapy_groups_without_ids(self) -> None:
        """Assign deterministic local IDs to referenceable group value objects."""
        record = Mock()
        record.model_dump.return_value = {
            "id": "civic.aid:1",
            "type": "Statement",
            "condition": {
                "membershipOperator": "AND",
                "conditions": [
                    {
                        "id": "civic.did:1",
                        "conceptType": "Disease",
                        "name": "Disease",
                    },
                    {
                        "id": "civic.phenotype:2",
                        "conceptType": "Phenotype",
                        "name": "Phenotype",
                    },
                ],
            },
            "therapeutic": {
                "membershipOperator": "AND",
                "therapies": [
                    {
                        "id": "civic.tid:1",
                        "conceptType": "Therapy",
                        "name": "Therapy 1",
                    },
                    {
                        "id": "civic.tid:2",
                        "conceptType": "Therapy",
                        "name": "Therapy 2",
                    },
                ],
            },
        }

        bundle = build_gks_bundle(
            [record],
            GksOutputMetadata(va_spec_python_version="test", created_at="2026-08-03"),
            [],
        )

        assert len(bundle.conditionSet) == 1
        assert len(bundle.therapyGroup) == 1
        condition_set_id = next(iter(bundle.conditionSet))
        therapy_group_id = next(iter(bundle.therapyGroup))
        assert condition_set_id.startswith("civic.conditionSet:")
        assert therapy_group_id.startswith("civic.therapyGroup:")
        assert bundle.conditionSet[condition_set_id]["id"] == condition_set_id
        assert bundle.therapyGroup[therapy_group_id]["id"] == therapy_group_id
        assert set(bundle.disease) == {"civic.did:1"}
        assert set(bundle.phenotype) == {"civic.phenotype:2"}
        assert bundle.conditionSet[condition_set_id]["conditions"] == [
            "#/disease/civic.did:1",
            "#/phenotype/civic.phenotype:2",
        ]
        statement = bundle.assertion["civic.aid:1"]
        assert statement["condition"] == f"#/conditionSet/{condition_set_id}"
        assert statement["therapeutic"] == f"#/therapyGroup/{therapy_group_id}"

    def test_rejects_invalid_group_without_id(self) -> None:
        """Reject malformed group-shaped data instead of leaving it inline."""
        record = Mock()
        record.model_dump.return_value = {
            "id": "civic.aid:1",
            "type": "Statement",
            "therapeutic": {"therapies": []},
        }

        with pytest.raises(ValidationError):
            build_gks_bundle(
                [record],
                GksOutputMetadata(
                    va_spec_python_version="test", created_at="2026-08-03"
                ),
                [],
            )
