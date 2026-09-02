"""Build referenced CIViC GKS bundles from inlined GKS Statements.

This module owns object routing, deduplication, reference creation, and computed
identifier assignment. The bundle's Pydantic schema remains in
:mod:`civicpy.exports.gks.bundle.models`.
"""

import json
import logging
from collections.abc import Iterable, Mapping
from datetime import date
from typing import Any

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.va_spec.base import (
    ConditionSet,
    CoreType,
    TherapyGroup,
    VariantClinicalSignificanceProposition,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from ga4gh.vrs.models import SequenceLocation, Syntax, VrsType
from civicpy.exports.gks.constants import (
    ALLELE_ORIGIN_QUALIFIER_FIELD,
    CODE_FIELD,
    CODING_FIELD,
    CONDITIONS_FIELD,
    GA4GH_CURIE_PREFIX,
    ID_FIELD,
    MAPPINGS_FIELD,
    PROPOSITION_FIELD,
    REFGET_ACCESSION_FIELD,
    TARGET_PROPOSITION_FIELD,
    THERAPIES_FIELD,
    TYPE_FIELD,
    CuriePrefix,
)
from civicpy.exports.gks.identifiers import (
    Identifiable,
    compute_identifier,
)
from civicpy.exports.gks.models import (
    GksAssertionError,
    GksOutputMetadata,
    GksRecord,
)
from civicpy.exports.gks.bundle.models import (
    Collection,
    GksBundle,
    GksBundleObject,
    GksBundleReference,
    Metadata,
    Statistics,
    VariantRepresentation,
)

_logger = logging.getLogger(__name__)


class BundleError(ValueError):
    """Indicate that an object cannot be represented in a GKS bundle."""


_PROPOSITION_FIELDS = frozenset({PROPOSITION_FIELD, TARGET_PROPOSITION_FIELD})
_VRS_VARIATION_TYPES = frozenset(
    {VrsType.ALLELE.value, VrsType.CN_CHANGE.value, VrsType.CN_COUNT.value}
)
_VRS_COPY_NUMBER_TYPES = frozenset(
    {VrsType.CN_CHANGE.value, VrsType.CN_COUNT.value}
)
_VRS_REPRESENTATION_BY_SYNTAX = {
    Syntax.HGVS_P.value: VariantRepresentation.PROTEIN,
    Syntax.HGVS_C.value: VariantRepresentation.CODING,
    Syntax.HGVS_G.value: VariantRepresentation.GENOMIC,
}
_GROUP_COLLECTION_BY_MEMBER_FIELD = {
    CONDITIONS_FIELD: Collection.CONDITION_SET,
    THERAPIES_FIELD: Collection.THERAPY_GROUP,
}

_BUNDLE_COLLECTION_BY_ID_PREFIX = {
    # SequenceReference uses a refget accession and has no VRS ``ga4gh.prefix``.
    f"{GA4GH_CURIE_PREFIX}:SQ.": Collection.SEQUENCE_REFERENCE,
    f"{GA4GH_CURIE_PREFIX}:{SequenceLocation.ga4gh.prefix}.": Collection.LOCATION,
    f"{CuriePrefix.GENE}:": Collection.FEATURE,
    f"{CuriePrefix.MOLECULAR_PROFILE}:": Collection.MOLECULAR_PROFILE,
    f"{CuriePrefix.DISEASE}:": Collection.DISEASE,
    f"{CuriePrefix.PHENOTYPE}:": Collection.PHENOTYPE,
    f"{CuriePrefix.THERAPY}:": Collection.THERAPY,
    f"{CuriePrefix.EVIDENCE}:": Collection.EVIDENCE,
    f"{CuriePrefix.ASSERTION}:": Collection.ASSERTION,
    f"{CuriePrefix.METHOD}:": Collection.METHOD,
    f"{CuriePrefix.ORGANIZATION}:": Collection.ORGANIZATION,
}

_BUNDLE_COLLECTION_BY_GKS_TYPE = {
    VrsType.SEQ_REF.value: Collection.SEQUENCE_REFERENCE,
    VrsType.SEQ_LOC.value: Collection.LOCATION,
    CategoricalVariant.model_fields[TYPE_FIELD].default: Collection.MOLECULAR_PROFILE,
    CoreType.DOCUMENT.value: Collection.SOURCE,
    CoreType.METHOD.value: Collection.METHOD,
    CoreType.AGENT.value: Collection.ORGANIZATION,
}

_COMPUTED_IDENTIFIER_MODELS: tuple[type[Identifiable], ...] = (
    VariantClinicalSignificanceProposition,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
_COMPUTED_IDENTIFIER_MODEL_BY_TYPE = {
    model.model_fields[TYPE_FIELD].default: model
    for model in _COMPUTED_IDENTIFIER_MODELS
}

_COMPUTED_IDENTIFIER_MODEL_BY_COLLECTION = {
    Collection.CONDITION_SET: ConditionSet,
    Collection.THERAPY_GROUP: TherapyGroup,
}


class _BundleBuilder:
    """Build a referenced bundle by extracting reusable nested GKS objects.

    Input Statements and output collection keys are sorted by identifier for
    deterministic serialization. If an identifier has multiple representations,
    the first is retained and a warning is logged.
    """

    def __init__(self) -> None:
        """Initialize empty collections and the deduplication index."""
        # Objects grouped by bundle collection.
        self._collections: dict[Collection, dict[str, GksBundleObject]] = {
            collection: {} for collection in Collection
        }
        # Original serialization for each collection key.
        self._source_serializations: dict[tuple[Collection, str], str] = {}
        # A molecular profile's CIViC mapping supplies the key for its defining VRS object.
        self._variant_reference_by_vrs_id: dict[str, GksBundleReference] = {}

    def build(
        self,
        records: Iterable[GksRecord],
        metadata: GksOutputMetadata,
        errors: list[GksAssertionError],
    ) -> GksBundle:
        """Build a validated bundle from inlined GKS Statements.

        :param records: Inlined VA-Spec GKS Statements.
        :param metadata: Provenance metadata for the generated bundle.
        :param errors: Transformation errors to include in the bundle.
        :return: Validated CIViC GKS bundle.
        """
        serialized_records = sorted(
            (record.model_dump(exclude_none=True) for record in records),
            key=lambda record: str(record.get(ID_FIELD, "")),
        )

        for record in serialized_records:
            collection = self._resolve_collection(record)
            if collection not in {
                Collection.EVIDENCE,
                Collection.ASSERTION,
            }:
                raise BundleError(
                    f"Unsupported Statement identifier: {record.get(ID_FIELD)!r}."
                )
            self._store_bundle_object(record, collection)

        bundle_metadata = Metadata(
            **metadata.model_dump(),
            statistics=Statistics.summarize_bundle_collections(
                self._collections
            ),
        )
        return GksBundle(
            **{
                collection.value: dict(sorted(objects.items()))
                for collection, objects in self._collections.items()
            },
            metadata=bundle_metadata,
            failed_assertion_ids=[error.assertion_id for error in errors],
            errors=errors,
        )

    @staticmethod
    def _build_reference(
        collection: Collection, collection_key: str
    ) -> GksBundleReference:
        """Build a JSON Pointer to an object in a bundle collection.

        :param collection: Root collection containing the object.
        :param collection_key: Key of the referenced object.
        :return: Document-local JSON Pointer.
        """
        return f"#/{collection.value}/{collection_key}"

    def _resolve_collection(
        self,
        record: Mapping[str, Any],
    ) -> Collection | None:
        """Select the root collection for a referenceable GKS object.

        Identifier namespaces take precedence over GKS types so semantic CIViC
        concepts remain in distinct collections.

        :param record: Serialized GKS object.
        :return: Matching collection, if the object can be referenced.
        """
        identifier = record.get(ID_FIELD)
        if isinstance(identifier, str):
            for prefix, collection in _BUNDLE_COLLECTION_BY_ID_PREFIX.items():
                if identifier.startswith(prefix):
                    return collection

        record_type = record.get(TYPE_FIELD)
        collection = (
            _BUNDLE_COLLECTION_BY_GKS_TYPE.get(record_type)
            if isinstance(record_type, str)
            else None
        )

        if collection is Collection.SEQUENCE_REFERENCE:
            refget_accession = record.get(REFGET_ACCESSION_FIELD)
            return collection if isinstance(refget_accession, str) else None

        return collection if isinstance(identifier, str) else None

    def _get_collection_key(
        self, record: Mapping[str, Any], collection: Collection
    ) -> str:
        """Return the key used to store an object in its bundle collection.

        :param record: Serialized GKS object.
        :param collection: Root collection where the object belongs.
        :return: Identifier used as the collection key.
        """
        identifier = record.get(ID_FIELD)
        if isinstance(identifier, str):
            return identifier

        if collection is Collection.SEQUENCE_REFERENCE:
            refget_accession = record.get(REFGET_ACCESSION_FIELD)
            if isinstance(refget_accession, str):
                return refget_accession

        raise BundleError(
            f"Object in {collection.value!r} has no usable collection key."
        )

    def _is_duplicate(
        self,
        collection: Collection,
        collection_key: str,
        source_value: Any,
    ) -> bool:
        """Record a source value, returning whether its key was already seen."""
        source_key = (collection, collection_key)
        serialization = json.dumps(
            source_value, sort_keys=True, separators=(",", ":"), default=str
        )
        existing_serialization = self._source_serializations.get(source_key)
        if existing_serialization is None:
            self._source_serializations[source_key] = serialization
            return False

        if existing_serialization != serialization:
            _logger.warning(
                "Multiple representations found for %s/%s; retaining the first "
                "bundle object.",
                collection.value,
                collection_key,
            )
        return True

    def _store_bundle_object(
        self, record: dict[str, Any], collection: Collection
    ) -> GksBundleReference:
        """Store a referenceable object once and return its JSON Pointer.

        Differing representations of a repeated identifier log a warning;
        the first representation is retained and later objects are not traversed.

        :param record: Serialized, identified GKS object to store.
        :param collection: Root collection where the object belongs.
        :return: JSON Pointer to the stored object.
        """
        # Build the collection key and a stable form used to detect conflicts.
        collection_key = self._get_collection_key(record, collection)
        if self._is_duplicate(collection, collection_key, record):
            return self._build_reference(collection, collection_key)

        if collection is Collection.MOLECULAR_PROFILE:
            self._store_molecular_profile_variants(record)

        # Replace nested objects with references, then store the result.
        transformed = {
            key: self._replace_nested_objects_with_references(value, key)
            for key, value in record.items()
        }
        self._collections[collection][collection_key] = transformed
        return self._build_reference(collection, collection_key)

    def _store_molecular_profile_variants(
        self, molecular_profile: Mapping[str, Any]
    ) -> None:
        """Group a profile's VRS representations under its CIViC VID."""
        variant_id = self._find_civic_variant_id(molecular_profile)
        if variant_id is None:
            return

        variations, defining_vrs_ids = self._collect_vrs_representations(
            molecular_profile
        )
        grouped = self._group_vrs_representations(
            variant_id, variations, defining_vrs_ids
        )
        if not grouped:
            original_variant = self._find_original_variant_coding(
                molecular_profile, variant_id
            )
            if original_variant is None:
                return

            unsupported = {
                key: self._replace_nested_objects_with_references(value, key)
                for key, value in original_variant.items()
            }
            self._store_variant_representations(
                variant_id,
                {VariantRepresentation.UNSUPPORTED.value: unsupported},
            )
            return

        transformed = {
            representation: {
                vrs_id: {
                    key: self._replace_nested_objects_with_references(value, key)
                    for key, value in variation.items()
                }
                for vrs_id, variation in representations.items()
            }
            for representation, representations in grouped.items()
        }

        self._store_variant_representations(variant_id, transformed)

    def _store_variant_representations(
        self, variant_id: str, representations: GksBundleObject
    ) -> None:
        """Store one CIViC variant's bundle representations.

        :param variant_id: CIViC variant collection key.
        :param representations: Variant representations to store.
        :return: ``None``.
        """
        if self._is_duplicate(Collection.VARIANT, variant_id, representations):
            return

        self._collections[Collection.VARIANT][variant_id] = representations

    def _collect_vrs_representations(
        self, molecular_profile: Mapping[str, Any]
    ) -> tuple[dict[str, dict[str, Any]], set[str]]:
        """Collect a profile's identified VRS objects and defining object IDs."""
        variations: dict[str, dict[str, Any]] = {}
        defining_vrs_ids: set[str] = set()
        for constraint in molecular_profile.get("constraints", []):
            if not isinstance(constraint, Mapping):
                continue
            for value in constraint.values():
                identifier = self._add_vrs_representation(variations, value)
                if identifier is not None:
                    defining_vrs_ids.add(identifier)

        for member in molecular_profile.get("members", []):
            self._add_vrs_representation(variations, member)

        return variations, defining_vrs_ids

    def _group_vrs_representations(
        self,
        variant_id: str,
        variations: Mapping[str, dict[str, Any]],
        defining_vrs_ids: set[str],
    ) -> dict[str, GksBundleObject]:
        """Group VRS objects by coordinate level and index their references."""
        grouped: dict[str, GksBundleObject] = {}
        for vrs_id, variation in variations.items():
            representation = self._classify_vrs_representation(
                variation, is_defining=vrs_id in defining_vrs_ids
            )
            reference = (
                f"#/{Collection.VARIANT.value}/{variant_id}/"
                f"{representation.value}/{vrs_id}"
            )
            self._variant_reference_by_vrs_id[vrs_id] = reference
            grouped.setdefault(representation.value, {})[vrs_id] = variation

        return grouped

    @staticmethod
    def _find_original_variant_coding(
        molecular_profile: Mapping[str, Any], variant_id: str
    ) -> dict[str, Any] | None:
        """Return the original variant coding mapped to a non-VRS profile.

        :param molecular_profile: Serialized Cat-VRS molecular profile.
        :param variant_id: CIViC variant identifier to find.
        :return: Original variant coding, if the profile maps to it.
        """
        for mapping in molecular_profile.get(MAPPINGS_FIELD, []):
            if not isinstance(mapping, Mapping):
                continue
            coding = mapping.get(CODING_FIELD)
            if not isinstance(coding, Mapping):
                continue
            if coding.get(ID_FIELD) == variant_id:
                return dict(coding)

        return None

    @staticmethod
    def _add_vrs_representation(
        variations: dict[str, dict[str, Any]], value: Any
    ) -> str | None:
        """Add an identified VRS variation and return its identifier."""
        if (
            not isinstance(value, Mapping)
            or value.get(TYPE_FIELD) not in _VRS_VARIATION_TYPES
        ):
            return None

        identifier = value.get(ID_FIELD)
        if isinstance(identifier, str):
            variations.setdefault(identifier, dict(value))
            return identifier

        return None

    @staticmethod
    def _classify_vrs_representation(
        variation: Mapping[str, Any], *, is_defining: bool
    ) -> VariantRepresentation:
        """Classify a VRS object by its HGVS syntax and constraint context."""
        syntaxes = {
            expression.get("syntax")
            for expression in variation.get("expressions", [])
            if isinstance(expression, Mapping)
        }
        for syntax, representation in _VRS_REPRESENTATION_BY_SYNTAX.items():
            if syntax in syntaxes:
                return representation

        if is_defining and variation.get(TYPE_FIELD) == VrsType.ALLELE.value:
            return VariantRepresentation.PROTEIN

        if variation.get(TYPE_FIELD) in _VRS_COPY_NUMBER_TYPES:
            return VariantRepresentation.GENOMIC

        return VariantRepresentation.OTHER

    @staticmethod
    def _find_civic_variant_id(molecular_profile: Mapping[str, Any]) -> str | None:
        """Return the single CIViC VID recorded in a profile's mappings.

        :param molecular_profile: Serialized Cat-VRS molecular profile.
        :raises BundleError: If mappings identify multiple CIViC variants.
        :return: CIViC VID, if one is mapped.
        """
        variant_ids: set[str] = set()
        for mapping in molecular_profile.get(MAPPINGS_FIELD, []):
            if not isinstance(mapping, Mapping):
                continue
            coding = mapping.get(CODING_FIELD)
            if not isinstance(coding, Mapping):
                continue
            identifier = coding.get(ID_FIELD)
            prefix = f"{CuriePrefix.VARIANT}:"
            if isinstance(identifier, str) and identifier.startswith(prefix):
                variant_ids.add(identifier)

        if len(variant_ids) > 1:
            raise BundleError(
                f"Molecular profile {molecular_profile.get(ID_FIELD)!r} maps to "
                f"multiple CIViC variants: {', '.join(sorted(variant_ids))}."
            )

        return next(iter(variant_ids), None)

    def _replace_nested_objects_with_references(
        self, value: Any, field_name: str
    ) -> Any:
        """Replace referenceable nested objects with JSON Pointers recursively.

        :param value: Nested value to process.
        :param field_name: Parent field name, used to identify propositions.
        :return: Referenced value, pointer, or unchanged scalar.
        """
        if isinstance(value, list):
            return [
                self._replace_nested_objects_with_references(item, field_name)
                for item in value
            ]

        if not isinstance(value, dict):
            return value

        identifier = value.get(ID_FIELD)
        if isinstance(identifier, str):
            variant_reference = self._variant_reference_by_vrs_id.get(identifier)
            if variant_reference is not None:
                return variant_reference

        if field_name in _PROPOSITION_FIELDS and ID_FIELD not in value:
            return self._store_object_with_computed_id(
                value, Collection.PROPOSITION
            )

        if field_name == ALLELE_ORIGIN_QUALIFIER_FIELD and ID_FIELD not in value:
            return self._store_variant_origin(value)

        collection = self._resolve_collection(value)
        if collection:
            return self._store_bundle_object(value, collection)

        group_collection = self._resolve_group_collection(value)
        if group_collection:
            return self._store_object_with_computed_id(value, group_collection)

        return {
            key: self._replace_nested_objects_with_references(item, key)
            for key, item in value.items()
        }

    @staticmethod
    def _resolve_group_collection(
        value: Mapping[str, Any],
    ) -> Collection | None:
        """Identify a group without a source-provided ID from its member field.

        :param value: Serialized nested GKS object.
        :return: Matching group collection, or ``None`` for another value object.
        """
        for member_field, collection in _GROUP_COLLECTION_BY_MEMBER_FIELD.items():
            if member_field in value:
                return collection

        return None

    def _store_object_with_computed_id(
        self,
        value: dict[str, Any],
        collection: Collection,
    ) -> GksBundleReference:
        """Store an object without a source ID under a computed CIViC GKS identifier.

        :param value: Serialized GKS object without a source-provided ID.
        :param collection: Root collection selected for the object.
        :return: JSON Pointer to the stored object.
        """
        gks_object = self._parse_computed_identifier_object(value, collection)
        identifier = compute_identifier(gks_object)
        return self._store_bundle_object({ID_FIELD: identifier, **value}, collection)

    def _store_variant_origin(self, value: dict[str, Any]) -> GksBundleReference:
        """Store a variant origin under the code from its first mapping."""
        try:
            code = value[MAPPINGS_FIELD][0][CODING_FIELD][CODE_FIELD]
        except (KeyError, IndexError, TypeError) as error:
            raise BundleError(
                "A variant origin requires a code in its first mapping."
            ) from error

        if not isinstance(code, str) or not code:
            raise BundleError(
                "A variant origin requires a code in its first mapping."
            )

        identifier = f"{CuriePrefix.VARIANT_ORIGIN}:{code}"
        return self._store_bundle_object(
            {ID_FIELD: identifier, **value}, Collection.VARIANT_ORIGIN
        )

    @staticmethod
    def _parse_computed_identifier_object(
        value: dict[str, Any],
        collection: Collection,
    ) -> Identifiable:
        """Parse an object that supports a computed CIViC GKS identifier.

        :param value: Serialized GKS object without a source-provided ID.
        :param collection: Bundle collection selected for the object.
        :raises BundleError: If the serialized object is unsupported.
        :return: Validated Pydantic GKS model.
        """
        model_class = _COMPUTED_IDENTIFIER_MODEL_BY_COLLECTION.get(collection)
        object_type = value.get(TYPE_FIELD)

        if model_class is None:
            model_class = (
                _COMPUTED_IDENTIFIER_MODEL_BY_TYPE.get(object_type)
                if isinstance(object_type, str)
                else None
            )

        if model_class is None:
            raise BundleError(
                f"Unsupported object for a computed CIViC GKS identifier: {object_type!r}."
            )

        model_value = value
        if TYPE_FIELD not in model_class.model_fields:
            model_value = {
                key: item for key, item in value.items() if key != TYPE_FIELD
            }

        return model_class.model_validate(model_value)


def build_gks_bundle(
    records: Iterable[GksRecord],
    metadata: GksOutputMetadata | None = None,
    errors: list[GksAssertionError] | None = None,
) -> GksBundle:
    """Build a referenced GKS bundle from inlined CIViC GKS Statements.

    :param records: Inlined VA-Spec GKS Statements translated from CIViC
        Assertions.
    :param metadata: Provenance metadata for the generated bundle. Defaults to
        metadata using today's date.
    :param errors: Transformation errors to include in the bundle. Defaults to
        an empty list.
    :return: Validated, reference-linked CIViC GKS bundle.
    """
    return _BundleBuilder().build(
        records,
        metadata or GksOutputMetadata(created_at=date.today().isoformat()),
        errors or [],
    )
