"""Compute stable CIViC identifiers for GKS objects without source IDs.

The CIViC GKS Bundle Format uses these identifiers for referenceable objects
without an ID defined by CIViC. For example, molecular profiles already have
``civic.mpid`` IDs, while propositions and groups need computed identifiers.
"""

import json
from typing import Any, TypeAlias

from ga4gh.core import sha512t24u
from ga4gh.core.models import ConceptMapping, MappableConcept
from ga4gh.va_spec.base import (
    ClinicalVariantProposition,
    ConditionSet,
    TherapyGroup,
)
from ga4gh.vrs.models import VrsType
from pydantic import BaseModel, Field

from civicpy.exports.gks.constants import (
    ALIASES_FIELD,
    CONDITIONS_FIELD,
    DESCRIPTION_FIELD,
    EXTENSIONS_FIELD,
    ID_FIELD,
    MAPPINGS_FIELD,
    MEMBERSHIP_OPERATOR_FIELD,
    NAME_FIELD,
    REFGET_ACCESSION_FIELD,
    THERAPIES_FIELD,
    TYPE_FIELD,
    CuriePrefix,
)


def _canonical_json(value: Any) -> str:
    """Serialize a JSON-compatible value deterministically."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


_COMMON_NON_IDENTITY_FIELDS = frozenset({ID_FIELD, EXTENSIONS_FIELD})
_NON_IDENTITY_PROPOSITION_FIELDS = _COMMON_NON_IDENTITY_FIELDS | frozenset(
    {NAME_FIELD, DESCRIPTION_FIELD, ALIASES_FIELD}
)


class IdentifierError(ValueError):
    """Indicate that an object does not support a computed CIViC GKS identifier."""


class AlleleOriginQualifier(MappableConcept):
    """Represent an allele origin concept with computable bundle identity.

    Both fields are required because they exclusively define the computed
    identifier for this GKS Bundle Format concept.
    """

    name: str
    mappings: list[ConceptMapping] = Field(
        min_length=1,
    )


# Object types supported by ``compute_identifier``.
Identifiable: TypeAlias = (
    ClinicalVariantProposition
    | ConditionSet
    | TherapyGroup
    | AlleleOriginQualifier
)


def compute_identifier(gks_object: BaseModel) -> str:
    """Compute a stable CIViC identifier for a supported GKS object.

    Supported objects are defined by ``Identifiable``.

    Identity is determined by these rules:

    * Existing ``id`` and ``extensions`` fields are excluded for all supported
      objects.
    * Group identity includes the membership operator and members. Identified
      members contribute only their IDs, and members are sorted because order is
      not semantically meaningful.
    * Proposition identity includes its type, predicate, and semantic participant
      or qualifier fields. ``name``, ``description``, and ``aliases`` are excluded.
      Nested condition sets and therapy groups contribute their own computed
      computed CIViC GKS identifiers.
    * Allele origin qualifier identity includes only ``name`` and ``mappings``.
      Mappings are sorted because their order is not semantically meaningful.
    * A nested ``SequenceReference`` contributes its ``refgetAccession`` when it
      has no ``id``.

    :param gks_object: Pydantic model to identify.
    :raises IdentifierError: If the object is not a supported GKS model.
    :return: Computed identifier in the namespace for the object's GKS type.
    """
    if isinstance(gks_object, ClinicalVariantProposition):
        identifier_prefix = CuriePrefix.PROPOSITION
        member_field = None
    elif isinstance(gks_object, ConditionSet):
        identifier_prefix = CuriePrefix.CONDITION_SET
        member_field = CONDITIONS_FIELD
    elif isinstance(gks_object, TherapyGroup):
        identifier_prefix = CuriePrefix.THERAPY_GROUP
        member_field = THERAPIES_FIELD
    elif isinstance(gks_object, AlleleOriginQualifier):
        identifier_prefix = CuriePrefix.VARIANT_ORIGIN
        member_field = MAPPINGS_FIELD
    else:
        raise IdentifierError(
            "Computed CIViC GKS identifiers are not supported for "
            f"{type(gks_object).__name__}."
        )

    if isinstance(gks_object, AlleleOriginQualifier):
        serialized_object = gks_object.model_dump(
            mode="json",
            include={NAME_FIELD, MAPPINGS_FIELD},
            exclude_none=True,
        )
    else:
        serialized_object = gks_object.model_dump(mode="json", exclude_none=True)
        excluded_fields = (
            _NON_IDENTITY_PROPOSITION_FIELDS
            if isinstance(gks_object, ClinicalVariantProposition)
            else _COMMON_NON_IDENTITY_FIELDS
        )
        serialized_object = {
            key: value
            for key, value in serialized_object.items()
            if key not in excluded_fields
        }
    identity_payload = _reduce_to_identity(
        serialized_object, use_computed_group_ids=False
    )

    if member_field and member_field in identity_payload:
        identity_payload[member_field] = sorted(
            identity_payload[member_field], key=_canonical_json
        )

    digest = sha512t24u(_canonical_json(identity_payload).encode())
    return f"{identifier_prefix}:{digest}"


def _reduce_to_identity(value: Any, *, use_computed_group_ids: bool = True) -> Any:
    """Reduce nested objects to the properties that define their identity.

    Identified nested objects reduce to ``id`` alone. A ``SequenceReference``
    without an ID reduces to ``refgetAccession``. Unidentified nested value
    objects retain their recursively projected semantic content.

    :param value: JSON-compatible value from a supported Pydantic GKS model.
    :param use_computed_group_ids: Whether nested groups without source IDs should
        reduce to their computed CIViC GKS identifiers. This is disabled for the root
        object.
    :return: Identity-defining content for deterministic hashing.
    """
    if isinstance(value, list):
        return [_reduce_to_identity(item) for item in value]

    if not isinstance(value, dict):
        return value

    identifier = value.get(ID_FIELD)
    if isinstance(identifier, str):
        return {ID_FIELD: identifier}

    refget_accession = value.get(REFGET_ACCESSION_FIELD)
    if value.get(TYPE_FIELD) == VrsType.SEQ_REF and isinstance(refget_accession, str):
        return {REFGET_ACCESSION_FIELD: refget_accession}

    if use_computed_group_ids:
        group_identifier = _compute_nested_group_identifier(value)

        if group_identifier:
            return {ID_FIELD: group_identifier}

    return {key: _reduce_to_identity(item) for key, item in value.items()}


def _compute_nested_group_identifier(value: dict[str, Any]) -> str | None:
    """Compute the identifier of a serialized nested GKS group when present.

    Upstream group models omit ``type`` when serialized, so required member and
    operator fields identify their concrete Pydantic class.

    :param value: Serialized nested value to inspect.
    :return: Computed CIViC GKS group identifier, or ``None`` for another
        value type.
    """
    if MEMBERSHIP_OPERATOR_FIELD not in value:
        return None

    if CONDITIONS_FIELD in value:
        group = ConditionSet.model_validate(value)
    elif THERAPIES_FIELD in value:
        group = TherapyGroup.model_validate(value)
    else:
        return None

    return compute_identifier(group)
