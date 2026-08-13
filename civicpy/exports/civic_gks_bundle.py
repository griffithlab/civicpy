"""Define the referenced CIViC GKS Bundle Format.

This module owns the public bundle schema, metadata, and statistics. Object
extraction and reference construction are implemented in
:mod:`civicpy.exports.civic_gks_bundle_builder`; file serialization remains in
:mod:`civicpy.exports.civic_gks_writer`.

The GKS Bundle Format permits pointers to every extracted collection. Organization
and proposition fields currently use pointers even though VA-Spec does not type
those fields as ``iriReference``.
"""

from collections import Counter
from collections.abc import Mapping
from enum import Enum
from types import MappingProxyType
from typing import Annotated, Any, TypeAlias

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.aac_2017 import (
    DiagnosticEvidenceLine,
    PrognosticEvidenceLine,
    TherapeuticEvidenceLine,
    VariantClinicalSignificanceStatement,
)
from ga4gh.va_spec.base import (
    Agent,
    ConditionSet,
    Document,
    Method,
    TherapyGroup,
    VariantClinicalSignificanceProposition,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from ga4gh.va_spec.ccv_2022 import (
    VariantOncogenicityEvidenceLine,
    VariantOncogenicityStatement,
)
from ga4gh.vrs.models import (
    Allele,
    CopyNumberChange,
    CopyNumberCount,
    SequenceLocation,
    SequenceReference,
)
from pydantic import BaseModel, Field, NonNegativeInt, SkipValidation, StringConstraints

from civicpy.exports.civic_gks_constants import (
    TYPE_FIELD,
    CivicGksBundleFormat,
)
from civicpy.exports.civic_gks_output import (
    GksAssertionError,
    GksOutputMetadata,
)

GksBundleObject: TypeAlias = dict[str, Any]
GksBundleReference: TypeAlias = str

# Bundle objects contain JSON Pointers in place of some nested model objects.
# SkipValidation preserves those referenced dictionaries while retaining the
# concrete upstream model in the generated JSON Schema.
GksVrsVariation: TypeAlias = Allele | CopyNumberChange | CopyNumberCount
GksProposition: TypeAlias = (
    VariantClinicalSignificanceProposition
    | VariantDiagnosticProposition
    | VariantOncogenicityProposition
    | VariantPrognosticProposition
    | VariantTherapeuticResponseProposition
)
GksEvidence: TypeAlias = (
    DiagnosticEvidenceLine
    | PrognosticEvidenceLine
    | TherapeuticEvidenceLine
    | VariantOncogenicityEvidenceLine
)
GksAssertion: TypeAlias = (
    VariantClinicalSignificanceStatement | VariantOncogenicityStatement
)
GksSequenceReferenceId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^SQ\.[A-Za-z0-9_-]+$")
]
GksLocationId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^ga4gh:SL\.[A-Za-z0-9_-]+$")
]
GksVariantId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.vid:[0-9]+$")
]
GksFeatureId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.gid:[0-9]+$")
]
GksMolecularProfileId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.mpid:[0-9]+$")
]
GksDiseaseId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.did:[0-9]+$")
]
GksPhenotypeId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.phenotype:[0-9]+$")
]
GksConditionSetId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.conditionSet:[A-Za-z0-9_-]+$")
]
GksTherapyId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.tid:[0-9]+$")
]
GksTherapyGroupId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.therapyGroup:[A-Za-z0-9_-]+$")
]
GksVariantOriginId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.variantOrigin:[A-Za-z0-9_-]+$")
]
GksSourceId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^(civic\.sid|pmid):[0-9]+$")
]
GksMethodId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.method:[A-Za-z0-9_-]+$")
]
GksOrganizationId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.organization:[A-Za-z0-9_-]+$")
]
GksPropositionId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.proposition:[A-Za-z0-9_-]+$")
]
GksEvidenceId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.eid:[0-9]+$")
]
GksAssertionId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.aid:[0-9]+$")
]
GksVariantRepresentations: TypeAlias = dict[
    str, dict[str, SkipValidation[GksVrsVariation]]
]

_CLOSED_COLLECTION_SCHEMA = {"additionalProperties": False}


class GksBundleCollection(str, Enum):
    """Names of keyed root collections in the CIViC GKS Bundle Format."""

    SEQUENCE_REFERENCE = "sequenceReference"
    LOCATION = "location"
    VARIANT = "variant"
    FEATURE = "feature"
    MOLECULAR_PROFILE = "molecularProfile"
    DISEASE = "disease"
    PHENOTYPE = "phenotype"
    CONDITION_SET = "conditionSet"
    THERAPY = "therapy"
    THERAPY_GROUP = "therapyGroup"
    VARIANT_ORIGIN = "variantOrigin"
    SOURCE = "source"
    METHOD = "method"
    ORGANIZATION = "organization"
    PROPOSITION = "proposition"
    EVIDENCE = "evidence"
    ASSERTION = "assertion"


class GksVariantRepresentation(str, Enum):
    """Coordinate levels used to group VRS representations of a CIViC variant."""

    PROTEIN = "protein"
    CODING = "coding"
    GENOMIC = "genomic"
    UNCLASSIFIED = "unclassified"


class CivicGksBundleError(ValueError):
    """Indicate that an object cannot be represented in a CIViC GKS bundle."""


# Include type counts only for collections that can contain more than one type.
_TYPE_COUNT_FIELD_BY_COLLECTION: Mapping[GksBundleCollection, str] = MappingProxyType(
    {
        GksBundleCollection.LOCATION: TYPE_FIELD,
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

            if collection is GksBundleCollection.VARIANT:
                type_counts = dict(
                    sorted(
                        Counter(
                            variation[TYPE_FIELD]
                            for representations in objects.values()
                            for variations in representations.values()
                            for variation in variations.values()
                        ).items()
                    )
                )
            elif type_field:
                type_counts = dict(
                    sorted(
                        Counter(
                            bundle_object[type_field]
                            for bundle_object in objects.values()
                        ).items()
                    )
                )
            else:
                type_counts = None

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
    become JSON Pointers such as ``#/molecularProfile/civic.mpid:33``.
    Propositions and groups without source-provided IDs receive deterministic
    bundle-local ``id`` values; other value objects remain inline.
    """

    sequenceReference: dict[
        GksSequenceReferenceId, SkipValidation[SequenceReference]
    ] = Field(json_schema_extra=_CLOSED_COLLECTION_SCHEMA)
    location: dict[GksLocationId, SkipValidation[SequenceLocation]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    variant: dict[GksVariantId, GksVariantRepresentations] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    feature: dict[GksFeatureId, SkipValidation[MappableConcept]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    molecularProfile: dict[
        GksMolecularProfileId, SkipValidation[CategoricalVariant]
    ] = Field(json_schema_extra=_CLOSED_COLLECTION_SCHEMA)
    disease: dict[GksDiseaseId, SkipValidation[MappableConcept]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    phenotype: dict[GksPhenotypeId, SkipValidation[MappableConcept]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    conditionSet: dict[GksConditionSetId, SkipValidation[ConditionSet]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    therapy: dict[GksTherapyId, SkipValidation[MappableConcept]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    therapyGroup: dict[GksTherapyGroupId, SkipValidation[TherapyGroup]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    variantOrigin: dict[GksVariantOriginId, SkipValidation[MappableConcept]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    source: dict[GksSourceId, SkipValidation[Document]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    method: dict[GksMethodId, SkipValidation[Method]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    organization: dict[GksOrganizationId, SkipValidation[Agent]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    proposition: dict[GksPropositionId, SkipValidation[GksProposition]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    evidence: dict[GksEvidenceId, SkipValidation[GksEvidence]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    assertion: dict[GksAssertionId, SkipValidation[GksAssertion]] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    metadata: GksBundleMetadata
    failed_assertion_ids: list[int]
    errors: list[GksAssertionError]
