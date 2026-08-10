"""Tests for reusable CIViC GKS computed identifiers."""

import pytest
from ga4gh.core.models import (
    Coding,
    ConceptMapping,
    MappableConcept,
    MembershipOperator,
    Relation,
    iriReference,
)
from ga4gh.va_spec.base import (
    ConditionSet,
    TherapyGroup,
    VariantClinicalSignificanceProposition,
    VariantOncogenicityProposition,
)
from pydantic import BaseModel

from civicpy.exports.civic_gks_identifier import (
    CivicGksAlleleOriginQualifier,
    CivicGksIdentifierError,
    compute_civic_gks_identifier,
)


def _allele_origin_qualifier(
    name: str, codes: list[str], concept_type: str | None = None
) -> CivicGksAlleleOriginQualifier:
    """Create an allele origin qualifier for computed-identifier tests.

    :param name: Allele origin display name included in identity.
    :param codes: Mapping codes included in identity.
    :param concept_type: Non-identity concept type used to test exclusion.
    :return: Allele origin qualifier with CIViC mappings.
    """
    return CivicGksAlleleOriginQualifier(
        name=name,
        conceptType=concept_type,
        mappings=[
            ConceptMapping(
                coding=Coding.model_validate(
                    {"code": code, "system": "https://civicdb.org"}
                ),
                relation=Relation.EXACT_MATCH,
            )
            for code in codes
        ],
    )


def test_allele_origin_identifier_uses_only_name_and_mappings() -> None:
    """Use name and mappings, but not other fields, for qualifier identity."""
    original = _allele_origin_qualifier("somatic", ["SOMATIC", "OTHER"])
    reordered = _allele_origin_qualifier(
        "somatic", ["OTHER", "SOMATIC"], concept_type="AlleleOrigin"
    ).model_copy(
        update={
            "id": "descriptive-system-id",
            "primaryCoding": Coding.model_validate(
                {"code": "IGNORED", "system": "https://example.org"}
            ),
        }
    )

    original_id = compute_civic_gks_identifier(original)

    assert original_id == compute_civic_gks_identifier(reordered)
    assert original_id != compute_civic_gks_identifier(
        _allele_origin_qualifier("germline", ["SOMATIC", "OTHER"])
    )
    assert original_id != compute_civic_gks_identifier(
        _allele_origin_qualifier("somatic", ["GERMLINE", "OTHER"])
    )


class _UnsupportedGksObject(BaseModel):
    """Represent an unsupported Pydantic model for error-path testing."""


def _concept(identifier: str, name: str) -> MappableConcept:
    """Create a minimally identified concept with descriptive content.

    :param identifier: Stable concept identifier.
    :param name: Descriptive concept name that must not define group identity.
    :return: Minimal GKS mappable concept.
    """
    return MappableConcept(id=identifier, conceptType="Therapy", name=name)


def test_group_identifier_uses_member_ids_not_descriptive_content() -> None:
    """Ignore member order and descriptive changes when identifying a group."""
    original = TherapyGroup(
        membershipOperator=MembershipOperator.AND,
        therapies=[
            _concept("civic.tid:1", "Therapy 1"),
            _concept("civic.tid:2", "Therapy 2"),
        ],
    )
    updated = TherapyGroup(
        id="civic.therapyGroup:preexisting",
        membershipOperator=MembershipOperator.AND,
        therapies=[
            _concept("civic.tid:2", "Renamed therapy 2"),
            _concept("civic.tid:1", "Renamed therapy 1"),
        ],
    )

    original_id = compute_civic_gks_identifier(original)

    assert original_id.startswith("civic.therapyGroup:")
    assert original_id == compute_civic_gks_identifier(updated)


def test_group_identifier_includes_membership_operator() -> None:
    """Treat different membership semantics as different group identities."""
    conditions: list[MappableConcept | ConditionSet] = [
        MappableConcept(id="civic.did:1", conceptType="Disease", name="Disease"),
        MappableConcept(
            id="civic.phenotype:2", conceptType="Phenotype", name="Phenotype"
        ),
    ]
    and_group = ConditionSet(
        membershipOperator=MembershipOperator.AND,
        conditions=conditions,
    )
    or_group = ConditionSet(
        membershipOperator=MembershipOperator.OR,
        conditions=conditions,
    )

    and_identifier = compute_civic_gks_identifier(and_group)

    assert and_identifier.startswith("civic.conditionSet:")
    assert and_identifier != compute_civic_gks_identifier(or_group)


def test_computes_proposition_identifier_from_pydantic_type() -> None:
    """Use proposition semantics while ignoring its descriptive metadata."""
    proposition = VariantOncogenicityProposition(
        name="Original name",
        description="Original description",
        subjectVariant=iriReference(root="civic.mpid:1"),
        objectTumorType=iriReference(root="civic.did:1"),
    )
    updated_proposition = proposition.model_copy(
        update={"name": "Updated name", "description": "Updated description"}
    )
    different_object = proposition.model_copy(
        update={"objectTumorType": iriReference(root="civic.did:2")}
    )

    identifier = compute_civic_gks_identifier(proposition)

    assert identifier.startswith("civic.proposition:")
    assert identifier == compute_civic_gks_identifier(updated_proposition)
    assert identifier != compute_civic_gks_identifier(different_object)


def test_proposition_identifier_uses_nested_group_identifier() -> None:
    """Ignore nested condition descriptions and order in proposition identity."""
    original_conditions = ConditionSet(
        membershipOperator=MembershipOperator.AND,
        conditions=[
            MappableConcept(id="civic.did:1", conceptType="Disease", name="Disease"),
            MappableConcept(
                id="civic.phenotype:2", conceptType="Phenotype", name="Phenotype"
            ),
        ],
    )
    updated_conditions = ConditionSet(
        membershipOperator=MembershipOperator.AND,
        conditions=[
            MappableConcept(
                id="civic.phenotype:2",
                conceptType="Phenotype",
                name="Renamed phenotype",
            ),
            MappableConcept(
                id="civic.did:1", conceptType="Disease", name="Renamed disease"
            ),
        ],
    )
    original = VariantClinicalSignificanceProposition.model_validate(
        {
            "subjectVariant": iriReference(root="civic.mpid:1"),
            "objectCondition": original_conditions,
        }
    )
    updated = VariantClinicalSignificanceProposition.model_validate(
        {
            "subjectVariant": iriReference(root="civic.mpid:1"),
            "objectCondition": updated_conditions,
        }
    )

    assert compute_civic_gks_identifier(original) == compute_civic_gks_identifier(
        updated
    )


def test_rejects_unsupported_pydantic_model() -> None:
    """Reject Pydantic objects whose identities are supplied by another model."""
    with pytest.raises(CivicGksIdentifierError, match="not supported"):
        compute_civic_gks_identifier(_UnsupportedGksObject())
