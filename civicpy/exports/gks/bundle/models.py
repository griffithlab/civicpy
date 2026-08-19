"""Define the referenced CIViC GKS Bundle Format.

This module owns the public bundle schema, metadata, and statistics. Object
extraction and reference construction are implemented in
:mod:`civicpy.exports.gks.bundle.builder`; file serialization remains in
:mod:`civicpy.exports.civic_gks_writer`.

The GKS Bundle Format permits pointers to every extracted collection. Organization
and proposition fields currently use pointers even though VA-Spec does not type
those fields as ``iriReference``.
"""

from collections import Counter
from collections.abc import Mapping
from enum import Enum
from typing import Annotated, Any, TypeAlias

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.aac_2017 import VariantClinicalSignificanceStatement
from ga4gh.va_spec.base import (
    Agent,
    ConditionSet,
    Document,
    Method,
    Statement,
    TherapyGroup,
    VariantClinicalSignificanceProposition,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
from ga4gh.vrs.models import (
    Allele,
    CopyNumberChange,
    SequenceLocation,
    SequenceReference,
)
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    StringConstraints,
    WithJsonSchema,
)

from civicpy.exports.gks.constants import (
    BUNDLE_FORMAT_NAME,
    BUNDLE_FORMAT_VERSION,
    TYPE_FIELD,
)
from civicpy.exports.gks.models import (
    GksAssertionError,
    GksModel,
    GksOutputMetadata,
)

GksBundleObject: TypeAlias = dict[str, Any]
GksBundleReference: TypeAlias = str
BUNDLE_SCHEMA_FILENAME = (
    f"{BUNDLE_FORMAT_NAME}-v{BUNDLE_FORMAT_VERSION}.schema.json"
)
BUNDLE_SCHEMA_ID = (
    f"urn:civic:gks-bundle:schema:{BUNDLE_FORMAT_VERSION}"
)


def _external_gks_schema(
    *models: type[BaseModel], concept_type: str | None = None
) -> WithJsonSchema:
    """Describe a value using canonical upstream GKS schema references."""
    references = [{"$ref": model.schema_id()} for model in models]
    schema: dict[str, Any] = (
        references[0] if len(references) == 1 else {"anyOf": references}
    )
    if concept_type is not None:
        schema = {
            "allOf": [
                schema,
                {
                    "properties": {"conceptType": {"const": concept_type}},
                    "required": ["conceptType"],
                },
            ]
        }
    return WithJsonSchema(schema)


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
_CLOSED_COLLECTION_SCHEMA = {"additionalProperties": False}


class Collection(str, Enum):
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


class VariantRepresentation(str, Enum):
    """Coordinate levels used to group VRS representations of a CIViC variant."""

    PROTEIN = "protein"
    CODING = "coding"
    GENOMIC = "genomic"
    UNCLASSIFIED = "unclassified"


def _variant_representations_schema() -> WithJsonSchema:
    """Describe the named coordinate levels containing VRS representations."""
    representation_collection = {
        "type": "object",
        "additionalProperties": {
            "anyOf": [
                {"$ref": Allele.schema_id()},
                {"$ref": CopyNumberChange.schema_id()},
            ]
        },
    }
    return WithJsonSchema(
        {
            "type": "object",
            "properties": {
                representation.value: representation_collection
                for representation in VariantRepresentation
            },
            "additionalProperties": False,
        }
    )


VariantRepresentations: TypeAlias = Annotated[
    dict[
        str,
        dict[
            str,
            Annotated[
                GksBundleObject,
                _external_gks_schema(Allele, CopyNumberChange),
            ],
        ],
    ],
    _variant_representations_schema(),
]


# Include type counts only for collections that can contain more than one type.
_TYPE_COUNT_FIELD_BY_COLLECTION = {
    Collection.LOCATION: TYPE_FIELD,
    Collection.PROPOSITION: TYPE_FIELD,
}


class CollectionStatistics(GksModel):
    """Summarize the objects stored in one root bundle collection."""

    count: NonNegativeInt = Field(description="Total objects in the collection.")
    types: dict[str, NonNegativeInt] | None = Field(
        default=None,
        description="Object counts grouped by concrete GKS type, when applicable.",
    )


class Statistics(GksModel):
    """Summarize the contents of a CIViC GKS bundle."""

    collections: dict[str, CollectionStatistics] = Field(
        description="Statistics keyed by root bundle collection name.",
    )

    @classmethod
    def summarize_bundle_collections(
        cls,
        collections: Mapping[Collection, Mapping[str, GksBundleObject]],
    ) -> "Statistics":
        """Calculate statistics from the bundle's root collections.

        :param collections: Bundle objects grouped by root collection.
        :return: Collection totals and applicable concrete-type breakdowns.
        """
        statistics: dict[str, CollectionStatistics] = {}
        for collection, objects in collections.items():
            type_field = _TYPE_COUNT_FIELD_BY_COLLECTION.get(collection)

            if collection is Collection.VARIANT:
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

            statistics[collection.value] = CollectionStatistics(
                count=len(objects),
                types=type_counts,
            )

        return cls(collections=statistics)


class Metadata(GksOutputMetadata):
    """Describe the bundle format, export provenance, and contents."""

    bundle_format: str = BUNDLE_FORMAT_NAME
    bundle_format_version: str = BUNDLE_FORMAT_VERSION
    statistics: Statistics = Field(
        description="Summary statistics for the bundle's root collections.",
    )


class GksBundle(GksModel):
    """Represent CIViC data in the GKS Bundle Format.

    Referenceable objects are stored once in keyed root collections. Nested uses
    become JSON Pointers such as ``#/molecularProfile/civic.mpid:33``.
    Propositions and groups without source-provided IDs receive deterministic
    bundle-local ``id`` values; other value objects remain inline.
    """

    model_config = ConfigDict(
        title=(f"CIViC GKS Bundle v{BUNDLE_FORMAT_VERSION}"),
        json_schema_extra={
            "$id": BUNDLE_SCHEMA_ID,
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "description": "CIViC data in the referenced GKS Bundle Format.",
            "civicBundleFormat": BUNDLE_FORMAT_NAME,
            "civicBundleFormatVersion": BUNDLE_FORMAT_VERSION,
        },
    )

    sequenceReference: dict[
        GksSequenceReferenceId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(SequenceReference),
        ],
    ] = Field(json_schema_extra=_CLOSED_COLLECTION_SCHEMA)
    location: dict[
        GksLocationId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(SequenceLocation),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    variant: dict[GksVariantId, VariantRepresentations] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    feature: dict[
        GksFeatureId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(MappableConcept, concept_type="Gene"),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    molecularProfile: dict[
        GksMolecularProfileId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(CategoricalVariant),
        ],
    ] = Field(json_schema_extra=_CLOSED_COLLECTION_SCHEMA)
    disease: dict[
        GksDiseaseId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(MappableConcept, concept_type="Disease"),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    phenotype: dict[
        GksPhenotypeId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(MappableConcept, concept_type="Phenotype"),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    conditionSet: dict[
        GksConditionSetId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(ConditionSet),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    therapy: dict[
        GksTherapyId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(MappableConcept, concept_type="Therapy"),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    therapyGroup: dict[
        GksTherapyGroupId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(TherapyGroup),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    variantOrigin: dict[
        GksVariantOriginId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(MappableConcept),
        ],
    ] = Field(json_schema_extra=_CLOSED_COLLECTION_SCHEMA)
    source: dict[
        GksSourceId,
        Annotated[GksBundleObject, _external_gks_schema(Document)],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    method: dict[
        GksMethodId,
        Annotated[GksBundleObject, _external_gks_schema(Method)],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    organization: dict[
        GksOrganizationId,
        Annotated[GksBundleObject, _external_gks_schema(Agent)],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    proposition: dict[
        GksPropositionId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(
                VariantClinicalSignificanceProposition,
                VariantDiagnosticProposition,
                VariantOncogenicityProposition,
                VariantPrognosticProposition,
                VariantTherapeuticResponseProposition,
            ),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    evidence: dict[
        GksEvidenceId,
        Annotated[GksBundleObject, _external_gks_schema(Statement)],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    assertion: dict[
        GksAssertionId,
        Annotated[
            GksBundleObject,
            _external_gks_schema(
                VariantClinicalSignificanceStatement,
                VariantOncogenicityStatement,
            ),
        ],
    ] = Field(
        json_schema_extra=_CLOSED_COLLECTION_SCHEMA
    )
    metadata: Metadata
    failed_assertion_ids: list[int]
    errors: list[GksAssertionError]
