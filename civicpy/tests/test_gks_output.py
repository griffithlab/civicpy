"""Tests for dereferenced and referenced CIViC GKS output documents."""

from unittest.mock import Mock

import pytest
from deepdiff import DeepDiff
from pydantic import ValidationError

from civicpy.civic import Assertion
from civicpy.exports.civic_gks_bundle import GksBundleOutput
from civicpy.exports.civic_gks_bundle_builder import build_gks_bundle
from civicpy.exports.civic_gks_output import (
    GksAssertionError,
    GksOutputMetadata,
)
from civicpy.exports.civic_gks_record import create_gks_record_from_assertion


class TestCivicGksBundleOutput:
    """Test the optional, reference-linked GKS Bundle Format export."""

    def test_schema_describes_concrete_bundle_objects(self) -> None:
        """Expose upstream GKS models instead of arbitrary JSON objects."""
        schema = GksBundleOutput.model_json_schema()
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

    def test_builds_expected_referenced_bundle(
        self,
        aid6: Assertion,
        aid202: Assertion,
        gks_bundle_expected: dict[str, object],
        mocked_normalizer: Mock,
    ) -> None:
        """Match the complete bundle representation assembled in ``conftest.py``."""
        actual_records = [
            create_gks_record_from_assertion(
                assertion,
            )
            for assertion in (aid6, aid202)
        ]
        metadata = GksOutputMetadata(
            va_spec_python_version="test", created_at="2026-08-03"
        )

        diff = DeepDiff(
            build_gks_bundle(actual_records, metadata, []).model_dump(
                exclude_none=True, serialize_as_any=True
            ),
            gks_bundle_expected,
            ignore_order=True,
        )
        assert diff == {}, diff

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
