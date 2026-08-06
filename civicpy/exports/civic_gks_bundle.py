"""Define the referenced CIViC GKS Bundle Format.

This module owns the public bundle schema, metadata, and statistics. Object
extraction and reference construction are implemented in
:mod:`civicpy.exports.civic_gks_bundle_builder`; file serialization remains in
:mod:`civicpy.exports.civic_gks_writer`.

The GKS Bundle Format permits pointers to every extracted collection. Agent and
proposition fields currently use pointers even though VA-Spec does not type
those fields as ``iriReference``.
"""

from collections import Counter
from collections.abc import Mapping
from enum import Enum
from types import MappingProxyType
from typing import Any, TypeAlias

from pydantic import BaseModel, Field, NonNegativeInt

from civicpy.exports.civic_gks_constants import (
    CONCEPT_TYPE_FIELD,
    TYPE_FIELD,
    CivicGksBundleFormat,
)
from civicpy.exports.civic_gks_output import (
    GksAssertionError,
    GksOutputMetadata,
)

GksBundleObject: TypeAlias = dict[str, Any]
GksBundleReference: TypeAlias = str


class GksBundleCollection(str, Enum):
    """Names of keyed root collections in the CIViC GKS Bundle Format."""

    SEQUENCE_REFERENCE = "sequenceReference"
    LOCATION = "location"
    MOLECULAR_VARIATION = "molecularVariation"
    GENE = "gene"
    CATEGORICAL_VARIANT = "categoricalVariant"
    CONDITION = "condition"
    CONDITION_SET = "conditionSet"
    THERAPY = "therapy"
    THERAPY_GROUP = "therapyGroup"
    ALLELE_ORIGIN_QUALIFIER = "alleleOriginQualifier"
    DOCUMENT = "document"
    METHOD = "method"
    AGENT = "agent"
    PROPOSITION = "proposition"
    STATEMENT = "statement"


class CivicGksBundleError(ValueError):
    """Indicate that an object cannot be represented in a CIViC GKS bundle."""


# Include type counts only for collections that can contain more than one type.
_TYPE_COUNT_FIELD_BY_COLLECTION: Mapping[GksBundleCollection, str] = MappingProxyType(
    {
        GksBundleCollection.LOCATION: TYPE_FIELD,
        GksBundleCollection.MOLECULAR_VARIATION: TYPE_FIELD,
        GksBundleCollection.CONDITION: CONCEPT_TYPE_FIELD,
        GksBundleCollection.PROPOSITION: TYPE_FIELD,
    }
)


class GksBundleCollectionStatistics(BaseModel):
    """Summarize the objects stored in one root bundle collection."""

    count: NonNegativeInt = Field(description="Total objects in the collection.")
    types: dict[str, NonNegativeInt] | None = Field(
        default=None,
        description="Object counts grouped by concrete GKS type, when applicable.",
    )


class GksBundleStatistics(BaseModel):
    """Summarize the contents of a CIViC GKS bundle."""

    collections: dict[str, GksBundleCollectionStatistics] = Field(
        description="Statistics keyed by root bundle collection name.",
    )

    @classmethod
    def summarize_bundle_collections(
        cls,
        collections: Mapping[GksBundleCollection, Mapping[str, GksBundleObject]],
    ) -> "GksBundleStatistics":
        """Calculate statistics from the bundle's root collections.

        :param collections: Bundle objects grouped by root collection.
        :return: Collection totals and applicable concrete-type breakdowns.
        """
        statistics: dict[str, GksBundleCollectionStatistics] = {}
        for collection, objects in collections.items():
            type_field = _TYPE_COUNT_FIELD_BY_COLLECTION.get(collection)

            type_counts = (
                dict(
                    sorted(
                        Counter(
                            bundle_object[type_field]
                            for bundle_object in objects.values()
                        ).items()
                    )
                )
                if type_field
                else None
            )

            statistics[collection.value] = GksBundleCollectionStatistics(
                count=len(objects),
                types=type_counts,
            )

        return cls(collections=statistics)


class GksBundleMetadata(GksOutputMetadata):
    """Describe the bundle format, export provenance, and contents."""

    bundle_format: str = CivicGksBundleFormat.NAME.value
    bundle_format_version: str = CivicGksBundleFormat.VERSION.value
    statistics: GksBundleStatistics = Field(
        description="Summary statistics for the bundle's root collections.",
    )


class GksBundleOutput(BaseModel):
    """Represent CIViC data in the GKS Bundle Format.

    Referenceable objects are stored once in keyed root collections. Nested uses
    become JSON Pointers such as ``#/categoricalVariant/civic.mpid:33``.
    Propositions and groups without source-provided IDs receive deterministic
    bundle-local ``id`` values; other value objects remain inline.
    """

    sequenceReference: dict[str, GksBundleObject]
    location: dict[str, GksBundleObject]
    molecularVariation: dict[str, GksBundleObject]
    gene: dict[str, GksBundleObject]
    categoricalVariant: dict[str, GksBundleObject]
    condition: dict[str, GksBundleObject]
    conditionSet: dict[str, GksBundleObject]
    therapy: dict[str, GksBundleObject]
    therapyGroup: dict[str, GksBundleObject]
    alleleOriginQualifier: dict[str, GksBundleObject]
    document: dict[str, GksBundleObject]
    method: dict[str, GksBundleObject]
    agent: dict[str, GksBundleObject]
    proposition: dict[str, GksBundleObject]
    statement: dict[str, GksBundleObject]
    metadata: GksBundleMetadata
    failed_assertion_ids: list[int]
    errors: list[GksAssertionError]
