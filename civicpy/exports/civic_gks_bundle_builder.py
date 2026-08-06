"""Build referenced CIViC GKS bundles from inlined GKS Statements.

This module owns object routing, deduplication, reference creation, and computed
identifier assignment. The bundle's Pydantic schema remains in
:mod:`civicpy.exports.civic_gks_bundle`.
"""

from collections.abc import Mapping
import logging
from types import MappingProxyType
from typing import Any, cast

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
from ga4gh.vrs.models import (
    Allele,
    CopyNumberChange,
    CopyNumberCount,
    SequenceLocation,
    VrsType,
)
from civicpy.exports.civic_gks_bundle import (
    CivicGksBundleError,
    GksBundleCollection,
    GksBundleMetadata,
    GksBundleObject,
    GksBundleOutput,
    GksBundleReference,
    GksBundleStatistics,
)
from civicpy.exports.civic_gks_constants import (
    ALLELE_ORIGIN_QUALIFIER_FIELD,
    CONDITIONS_FIELD,
    GA4GH_CURIE_PREFIX,
    ID_FIELD,
    PROPOSITION_FIELD,
    REFGET_ACCESSION_FIELD,
    TARGET_PROPOSITION_FIELD,
    THERAPIES_FIELD,
    TYPE_FIELD,
    CivicGksCuriePrefix,
)
from civicpy.exports.civic_gks_identifier import (
    CivicGksAlleleOriginQualifier,
    CivicGksComputedIdentifierObject,
    compute_civic_gks_identifier,
)
from civicpy.exports.civic_gks_output import (
    GksAssertionError,
    GksOutputMetadata,
    GksRecord,
)
from civicpy.exports.civic_gks_utils import serialize_canonical_json

_logger = logging.getLogger(__name__)

_PROPOSITION_FIELDS = frozenset({PROPOSITION_FIELD, TARGET_PROPOSITION_FIELD})

_BUNDLE_COLLECTION_BY_ID_PREFIX: Mapping[str, GksBundleCollection] = MappingProxyType(
    {
        # SequenceReference uses a refget accession and has no VRS ``ga4gh.prefix``.
        f"{GA4GH_CURIE_PREFIX}:SQ.": GksBundleCollection.SEQUENCE_REFERENCE,
        f"{GA4GH_CURIE_PREFIX}:{SequenceLocation.ga4gh.prefix}.": GksBundleCollection.LOCATION,
        f"{GA4GH_CURIE_PREFIX}:{Allele.ga4gh.prefix}.": GksBundleCollection.MOLECULAR_VARIATION,
        f"{GA4GH_CURIE_PREFIX}:{CopyNumberChange.ga4gh.prefix}.": GksBundleCollection.MOLECULAR_VARIATION,
        f"{GA4GH_CURIE_PREFIX}:{CopyNumberCount.ga4gh.prefix}.": GksBundleCollection.MOLECULAR_VARIATION,
        f"{CivicGksCuriePrefix.GENE.value}:": GksBundleCollection.GENE,
        f"{CivicGksCuriePrefix.MOLECULAR_PROFILE.value}:": GksBundleCollection.CATEGORICAL_VARIANT,
        f"{CivicGksCuriePrefix.DISEASE.value}:": GksBundleCollection.CONDITION,
        f"{CivicGksCuriePrefix.PHENOTYPE.value}:": GksBundleCollection.CONDITION,
        f"{CivicGksCuriePrefix.THERAPY.value}:": GksBundleCollection.THERAPY,
        f"{CivicGksCuriePrefix.EVIDENCE.value}:": GksBundleCollection.STATEMENT,
        f"{CivicGksCuriePrefix.ASSERTION.value}:": GksBundleCollection.STATEMENT,
        f"{CivicGksCuriePrefix.METHOD.value}:": GksBundleCollection.METHOD,
        f"{CivicGksCuriePrefix.ORGANIZATION.value}:": GksBundleCollection.AGENT,
    }
)

_BUNDLE_COLLECTION_BY_GKS_TYPE: Mapping[str, GksBundleCollection] = MappingProxyType(
    {
        VrsType.SEQ_REF.value: GksBundleCollection.SEQUENCE_REFERENCE,
        VrsType.SEQ_LOC.value: GksBundleCollection.LOCATION,
        VrsType.ALLELE.value: GksBundleCollection.MOLECULAR_VARIATION,
        VrsType.CN_CHANGE.value: GksBundleCollection.MOLECULAR_VARIATION,
        VrsType.CN_COUNT.value: GksBundleCollection.MOLECULAR_VARIATION,
        CategoricalVariant.model_fields[
            TYPE_FIELD
        ].default: GksBundleCollection.CATEGORICAL_VARIANT,
        CoreType.DOCUMENT.value: GksBundleCollection.DOCUMENT,
        CoreType.METHOD.value: GksBundleCollection.METHOD,
        CoreType.AGENT.value: GksBundleCollection.AGENT,
    }
)

_COMPUTED_IDENTIFIER_MODELS: tuple[type[CivicGksComputedIdentifierObject], ...] = (
    VariantClinicalSignificanceProposition,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
_COMPUTED_IDENTIFIER_MODEL_BY_TYPE: Mapping[
    str, type[CivicGksComputedIdentifierObject]
] = (
    MappingProxyType(
        {
            model.model_fields[TYPE_FIELD].default: model
            for model in _COMPUTED_IDENTIFIER_MODELS
        }
    )
)

_COMPUTED_IDENTIFIER_MODEL_BY_COLLECTION: Mapping[
    GksBundleCollection, type[CivicGksComputedIdentifierObject]
] = MappingProxyType(
    {
        GksBundleCollection.CONDITION_SET: ConditionSet,
        GksBundleCollection.THERAPY_GROUP: TherapyGroup,
        GksBundleCollection.ALLELE_ORIGIN_QUALIFIER: CivicGksAlleleOriginQualifier,
    }
)


class _BundleBuilder:
    """Build a referenced bundle by extracting reusable nested GKS objects.

    Input Statements and output collection keys are sorted by identifier for
    deterministic serialization. If an identifier has multiple representations,
    the first is retained and a warning is logged.
    """

    def __init__(self) -> None:
        """Initialize empty collections and the deduplication index."""
        # Objects grouped by bundle collection.
        self._collections: dict[GksBundleCollection, dict[str, GksBundleObject]] = {
            collection: {} for collection in GksBundleCollection
        }
        # Original serialization for each collection key.
        self._source_serializations: dict[tuple[GksBundleCollection, str], str] = {}

    def build(
        self,
        records: list[GksRecord],
        metadata: GksOutputMetadata,
        errors: list[GksAssertionError],
    ) -> GksBundleOutput:
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
            self._store_bundle_object(record, GksBundleCollection.STATEMENT)

        bundle_metadata = GksBundleMetadata(
            **metadata.model_dump(),
            statistics=GksBundleStatistics.summarize_bundle_collections(
                self._collections
            ),
        )
        return GksBundleOutput(
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
        collection: GksBundleCollection, collection_key: str
    ) -> GksBundleReference:
        """Build a JSON Pointer to an object in a bundle collection.

        :param collection: Root collection containing the object.
        :param collection_key: Key of the referenced object.
        :return: Document-local JSON Pointer.
        """
        return f"#/{collection.value}/{collection_key}"

    @staticmethod
    def _resolve_collection(
        record: Mapping[str, Any],
    ) -> GksBundleCollection | None:
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

        if collection is GksBundleCollection.SEQUENCE_REFERENCE:
            refget_accession = record.get(REFGET_ACCESSION_FIELD)
            return collection if isinstance(refget_accession, str) else None

        return collection if isinstance(identifier, str) else None

    @staticmethod
    def _get_collection_key(
        record: Mapping[str, Any], collection: GksBundleCollection
    ) -> str:
        """Return the key used to store an object in its bundle collection.

        :param record: Serialized GKS object.
        :param collection: Root collection where the object belongs.
        :return: Identifier used as the collection key.
        """
        identifier = record.get(ID_FIELD)
        if (
            not isinstance(identifier, str)
            and collection is GksBundleCollection.SEQUENCE_REFERENCE
        ):
            return cast(str, record[REFGET_ACCESSION_FIELD])

        return cast(str, identifier)

    def _store_bundle_object(
        self, record: dict[str, Any], collection: GksBundleCollection
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
        source_key = (collection, collection_key)
        stable_record = serialize_canonical_json(record)
        existing_source = self._source_serializations.get(source_key)

        # Reuse the first object found for an identifier.
        if existing_source is not None:
            if existing_source != stable_record:
                _logger.warning(
                    "Multiple representations found for %s/%s; retaining the first "
                    "bundle object.",
                    collection.value,
                    collection_key,
                )
            return self._build_reference(collection, collection_key)

        # Mark the object as seen before processing its nested values.
        self._source_serializations[source_key] = stable_record

        # Replace nested objects with references, then store the result.
        transformed = {
            key: self._replace_nested_objects_with_references(value, key)
            for key, value in record.items()
        }
        self._collections[collection][collection_key] = transformed
        return self._build_reference(collection, collection_key)

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

        if field_name in _PROPOSITION_FIELDS and ID_FIELD not in value:
            return self._store_object_with_computed_id(
                value, GksBundleCollection.PROPOSITION
            )

        if field_name == ALLELE_ORIGIN_QUALIFIER_FIELD and ID_FIELD not in value:
            return self._store_object_with_computed_id(
                value, GksBundleCollection.ALLELE_ORIGIN_QUALIFIER
            )

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
    ) -> GksBundleCollection | None:
        """Identify a group without a source-provided ID from its member field.

        :param value: Serialized nested GKS object.
        :return: Matching group collection, or ``None`` for another value object.
        """
        if CONDITIONS_FIELD in value:
            return GksBundleCollection.CONDITION_SET

        if THERAPIES_FIELD in value:
            return GksBundleCollection.THERAPY_GROUP

        return None

    def _store_object_with_computed_id(
        self,
        value: dict[str, Any],
        collection: GksBundleCollection,
    ) -> GksBundleReference:
        """Store an object without a source ID under a computed ``civic.gks`` identifier.

        :param value: Serialized GKS object without a source-provided ID.
        :param collection: Root collection selected for the object.
        :return: JSON Pointer to the stored object.
        """
        gks_object = self._parse_computed_identifier_object(value, collection)
        identifier = compute_civic_gks_identifier(gks_object)
        return self._store_bundle_object({ID_FIELD: identifier, **value}, collection)

    @staticmethod
    def _parse_computed_identifier_object(
        value: dict[str, Any],
        collection: GksBundleCollection,
    ) -> CivicGksComputedIdentifierObject:
        """Parse an object that supports a computed ``civic.gks`` identifier.

        :param value: Serialized GKS object without a source-provided ID.
        :param collection: Bundle collection selected for the object.
        :raises CivicGksBundleError: If the serialized object is unsupported.
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
            raise CivicGksBundleError(
                f"Unsupported object for a computed civic.gks identifier: {object_type!r}."
            )

        model_value = value
        if TYPE_FIELD not in model_class.model_fields:
            model_value = {
                key: item for key, item in value.items() if key != TYPE_FIELD
            }

        return model_class.model_validate(model_value)


def build_gks_bundle(
    records: list[GksRecord],
    metadata: GksOutputMetadata,
    errors: list[GksAssertionError],
) -> GksBundleOutput:
    """Build a referenced GKS bundle from inlined CIViC GKS Statements.

    :param records: Inlined VA-Spec GKS Statements translated from CIViC
        Assertions.
    :param metadata: Provenance metadata for the generated bundle.
    :param errors: Transformation errors to include in the bundle.
    :return: Validated, reference-linked CIViC GKS bundle.
    """
    return _BundleBuilder().build(records, metadata, errors)
