"""Define the referenced CIViC GKM Bundle Format.

This module owns the public bundle schema, metadata, and statistics. Object
extraction and reference construction are implemented in
:mod:`civicpy.exports.gks.bundle.builder`; file serialization remains in
:mod:`civicpy.exports.civic_gks_writer`.

The GKM Bundle Format permits pointers to every extracted collection. Organization
and proposition fields currently use pointers even though VA-Spec does not type
those fields as ``iriReference``.
"""

from collections import Counter
from collections.abc import Mapping
from enum import Enum
from typing import Annotated, Any, TypeAlias

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import Coding, MappableConcept
from ga4gh.va_spec.aac_2017 import VariantClinicalSignificanceStatement
from ga4gh.va_spec.base import (
    Agent,
    Condition,
    ConditionSet,
    Document,
    Method,
    Statement,
    Therapy,
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
    GkmAssertionError,
    GkmModel,
    GkmOutputMetadata,
)

GkmBundleObject: TypeAlias = dict[str, Any]
GkmBundleReference: TypeAlias = str
BUNDLE_SCHEMA_FILENAME = f"{BUNDLE_FORMAT_NAME}-v{BUNDLE_FORMAT_VERSION}.schema.json"
BUNDLE_SCHEMA_ID = f"urn:civic:gkm-bundle:schema:{BUNDLE_FORMAT_VERSION}"
CIVIC_KNOWLEDGE_MODEL_URL = "https://civic.readthedocs.io/en/latest/model.html"


def _external_gkm_schema(
    *models: type[BaseModel], concept_type: str | None = None
) -> WithJsonSchema:
    """Describe a value using canonical upstream GKM schema references."""
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


GkmSequenceReferenceId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^SQ\.[A-Za-z0-9_-]+$")
]
GkmLocationId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^ga4gh:SL\.[A-Za-z0-9_-]+$")
]
GkmVariantId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.vid:[0-9]+$")
]
GkmFeatureId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.gid:[0-9]+$")
]
GkmMolecularProfileId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.mpid:[0-9]+$")
]
GkmDiseaseId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.did:[0-9]+$")
]
GkmPhenotypeId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.phenotype:[0-9]+$")
]
GkmConditionSetId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.conditionSet:[A-Za-z0-9_-]+$")
]
GkmTherapyId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.tid:[0-9]+$")
]
GkmTherapyGroupId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.therapyGroup:[A-Za-z0-9_-]+$")
]
GkmVariantOriginId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.variantOrigin:[A-Za-z0-9_-]+$")
]
GkmSourceId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^(civic\.sid|pmid):[0-9]+$")
]
GkmMethodId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.method:[A-Za-z0-9_-]+$")
]
GkmOrganizationId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.organization:[A-Za-z0-9_-]+$")
]
GkmPropositionId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.proposition:[A-Za-z0-9_-]+$")
]
GkmEvidenceId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.eid:[0-9]+$")
]
GkmAssertionId: TypeAlias = Annotated[
    str, StringConstraints(pattern=r"^civic\.aid:[0-9]+$")
]


class Collection(str, Enum):
    """Names of keyed root collections in the CIViC GKM Bundle Format."""

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
    """Representation levels used to group CIViC variants."""

    PROTEIN = "protein"
    CODING = "coding"
    GENOMIC = "genomic"
    OTHER = "other"
    UNSUPPORTED = "unsupported"


_COLLECTION_DESCRIPTIONS: Mapping[Collection, str] = {
    Collection.SEQUENCE_REFERENCE: (
        "Reference sequences used by variant coordinates and HGVS expressions."
    ),
    Collection.LOCATION: (
        "Representative sequence locations for CIViC Variant coordinates."
    ),
    Collection.VARIANT: (
        "User-defined genomic alterations or related molecular events with "
        "potential clinical relevance."
    ),
    Collection.FEATURE: (
        "Genes associated with CIViC Variants and Molecular Profiles."
    ),
    Collection.MOLECULAR_PROFILE: (
        "CIViC Molecular Profiles interpreted for clinical relevance; this "
        "export currently includes profiles with one CIViC Variant."
    ),
    Collection.DISEASE: (
        "Cancer types or subtypes associated with CIViC Evidence Items and Assertions."
    ),
    Collection.PHENOTYPE: (
        "Symptoms or abnormalities from the Human Phenotype Ontology that add "
        "clinical context."
    ),
    Collection.CONDITION_SET: (
        "Disease and phenotype context grouped for reuse in GKM propositions."
    ),
    Collection.THERAPY: (
        "Drugs or treatment types associated with predictive Evidence Items and "
        "Assertions."
    ),
    Collection.THERAPY_GROUP: (
        "Groups of multiple therapies used in predictive interpretations."
    ),
    Collection.VARIANT_ORIGIN: (
        "Presumed source of a variant in the study context, such as somatic, "
        "germline, unknown, or not applicable."
    ),
    Collection.SOURCE: (
        "Publications or abstracts cited as support for CIViC Evidence Items "
        "and Assertions."
    ),
    Collection.METHOD: (
        "Curation or classification method used to produce a CIViC GKM Statement."
    ),
    Collection.ORGANIZATION: (
        "CIViC organizations connected to assertion approvals or submissions."
    ),
    Collection.PROPOSITION: (
        "Structured clinical interpretation evaluated by CIViC Evidence Items "
        "or Assertions."
    ),
    Collection.EVIDENCE: (
        "Clinical evidence statements manually curated from a source publication."
    ),
    Collection.ASSERTION: (
        "CIViC Assertions that summarize and classify evidence for a Molecular "
        "Profile in a disease context."
    ),
}
_COLLECTION_KEY_DESCRIPTIONS: Mapping[Collection, str] = {
    Collection.SEQUENCE_REFERENCE: "refget sequence accessions",
    Collection.LOCATION: "GA4GH SequenceLocation identifiers",
    Collection.VARIANT: "CIViC Variant IDs",
    Collection.FEATURE: "CIViC Gene IDs",
    Collection.MOLECULAR_PROFILE: "CIViC Molecular Profile IDs",
    Collection.DISEASE: "CIViC Disease IDs",
    Collection.PHENOTYPE: "CIViC Phenotype IDs",
    Collection.CONDITION_SET: "computed bundle-local condition set identifiers",
    Collection.THERAPY: "CIViC Therapy IDs",
    Collection.THERAPY_GROUP: "computed bundle-local therapy group identifiers",
    Collection.VARIANT_ORIGIN: "computed identifiers based on CIViC variant origin codes",
    Collection.SOURCE: "CIViC Source IDs or PubMed IDs",
    Collection.METHOD: "CIViC curation method IDs, such as civic.method:2019",
    Collection.ORGANIZATION: "CIViC Organization IDs",
    Collection.PROPOSITION: "computed bundle-local proposition identifiers",
    Collection.EVIDENCE: "CIViC Evidence IDs",
    Collection.ASSERTION: "CIViC Assertion IDs",
}
_VARIANT_REPRESENTATION_DESCRIPTIONS: Mapping[VariantRepresentation, str] = {
    VariantRepresentation.PROTEIN: (
        "HGVS protein expressions and other protein-level representations."
    ),
    VariantRepresentation.CODING: (
        "HGVS coding DNA expressions, usually relative to a transcript."
    ),
    VariantRepresentation.GENOMIC: (
        "HGVS genomic expressions and representative genomic coordinates."
    ),
    VariantRepresentation.OTHER: (
        "VRS representations that are not clearly protein, coding, or genomic."
    ),
    VariantRepresentation.UNSUPPORTED: (
        "Original variant identifier and label used when current tooling "
        "cannot represent the variant using VRS."
    ),
}


def _closed_collection_schema(collection: Collection) -> dict[str, Any]:
    """Return schema extras for a closed keyed bundle collection.

    :param collection: Bundle collection to describe.
    :return: JSON Schema extras that close the collection and describe its
        contents.
    """
    return {
        "additionalProperties": False,
        "description": (
            f"{_COLLECTION_DESCRIPTIONS[collection]} Collection keys are "
            f"{_COLLECTION_KEY_DESCRIPTIONS[collection]}."
        ),
    }


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
    properties = {
        representation.value: {
            **representation_collection,
            "description": _VARIANT_REPRESENTATION_DESCRIPTIONS[representation],
        }
        for representation in VariantRepresentation
        if representation is not VariantRepresentation.UNSUPPORTED
    }
    properties[VariantRepresentation.UNSUPPORTED.value] = {
        "$ref": Coding.schema_id(),
        "description": _VARIANT_REPRESENTATION_DESCRIPTIONS[
            VariantRepresentation.UNSUPPORTED
        ],
    }

    return WithJsonSchema(
        {
            "type": "object",
            "properties": properties,
            "additionalProperties": False,
        }
    )


VariantRepresentations: TypeAlias = Annotated[
    dict[str, Any],
    _variant_representations_schema(),
]


def _count_variant_representation_types(
    variants: Mapping[str, GkmBundleObject],
) -> dict[str, int]:
    """Count concrete object types stored in variant representations.

    :param variants: Variant representations keyed by CIViC variant ID.
    :return: Counts keyed by representation object type.
    """
    type_counts: Counter[str] = Counter()
    unsupported_field = VariantRepresentation.UNSUPPORTED.value

    for representations in variants.values():
        if unsupported_field in representations:
            type_counts["Coding"] += 1

        for representation, variations in representations.items():
            if representation == unsupported_field:
                continue

            for variation in variations.values():
                type_counts[variation[TYPE_FIELD]] += 1

    return dict(sorted(type_counts.items()))


# Include type counts only for collections that can contain more than one type.
_TYPE_COUNT_FIELD_BY_COLLECTION = {
    Collection.LOCATION: TYPE_FIELD,
    Collection.PROPOSITION: TYPE_FIELD,
}


class CollectionStatistics(GkmModel):
    """Summarize the objects stored in one root bundle collection."""

    count: NonNegativeInt = Field(description="Total objects in the collection.")
    types: dict[str, NonNegativeInt] | None = Field(
        default=None,
        description="Object counts grouped by concrete GKM type, when applicable.",
    )


class Statistics(GkmModel):
    """Summarize the contents of a CIViC GKM bundle."""

    collections: dict[str, CollectionStatistics] = Field(
        description="Statistics keyed by root bundle collection name.",
    )

    @classmethod
    def summarize_bundle_collections(
        cls,
        collections: Mapping[Collection, Mapping[str, GkmBundleObject]],
    ) -> "Statistics":
        """Calculate statistics from the bundle's root collections.

        :param collections: Bundle objects grouped by root collection.
        :return: Collection totals and applicable concrete-type breakdowns.
        """
        statistics: dict[str, CollectionStatistics] = {}
        for collection, objects in collections.items():
            type_field = _TYPE_COUNT_FIELD_BY_COLLECTION.get(collection)

            if collection is Collection.VARIANT:
                type_counts = _count_variant_representation_types(objects)
                count = sum(len(values) for values in objects.values())
            elif type_field:
                type_counts = dict(
                    sorted(
                        Counter(
                            bundle_object[type_field]
                            for bundle_object in objects.values()
                        ).items()
                    )
                )
                count = len(objects)
            else:
                type_counts = None
                count = len(objects)

            statistics[collection.value] = CollectionStatistics(
                count=count,
                types=type_counts,
            )

        return cls(collections=statistics)


class Metadata(GkmOutputMetadata):
    """Describe the bundle format, export provenance, and contents."""

    bundle_format: str = BUNDLE_FORMAT_NAME
    bundle_format_version: str = BUNDLE_FORMAT_VERSION
    statistics: Statistics = Field(
        description="Summary statistics for the bundle's root collections.",
    )


class GkmBundle(GkmModel):
    """Represent CIViC data in the GKM Bundle Format.

    Referenceable objects are stored once in keyed root collections. Nested uses
    become JSON Pointers such as ``#/molecularProfile/civic.mpid:33``.
    Propositions and groups without source-provided IDs receive deterministic
    bundle-local ``id`` values; other value objects remain inline.
    """

    model_config = ConfigDict(
        title=(f"CIViC GKM Bundle v{BUNDLE_FORMAT_VERSION}"),
        json_schema_extra={
            "$id": BUNDLE_SCHEMA_ID,
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "description": (
                "CIViC organizes curated cancer variant knowledge around "
                "variants and molecular profiles with evidence from source "
                "publications, summary assertions, and clinical context such "
                "as diseases, therapies, phenotypes, variant origins, and "
                "curating organizations. This schema describes a GA4GH GKM "
                "representation of those CIViC concepts. Top-level keys use "
                "CIViC knowledge model terms where possible, with additional "
                "keys for supporting GKM representation details. For the "
                f"CIViC data model, see {CIVIC_KNOWLEDGE_MODEL_URL}."
            ),
            "civicBundleFormat": BUNDLE_FORMAT_NAME,
            "civicBundleFormatVersion": BUNDLE_FORMAT_VERSION,
        },
    )

    sequenceReference: dict[
        GkmSequenceReferenceId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(SequenceReference),
        ],
    ] = Field(
        json_schema_extra=_closed_collection_schema(Collection.SEQUENCE_REFERENCE)
    )
    location: dict[
        GkmLocationId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(SequenceLocation),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.LOCATION))
    variant: dict[GkmVariantId, VariantRepresentations] = Field(
        json_schema_extra=_closed_collection_schema(Collection.VARIANT)
    )
    feature: dict[
        GkmFeatureId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(MappableConcept, concept_type="Gene"),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.FEATURE))
    molecularProfile: dict[
        GkmMolecularProfileId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(CategoricalVariant),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.MOLECULAR_PROFILE))
    disease: dict[
        GkmDiseaseId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(Condition, concept_type="Disease"),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.DISEASE))
    phenotype: dict[
        GkmPhenotypeId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(Condition, concept_type="Phenotype"),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.PHENOTYPE))
    conditionSet: dict[
        GkmConditionSetId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(ConditionSet),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.CONDITION_SET))
    therapy: dict[
        GkmTherapyId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(Therapy, concept_type="Therapy"),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.THERAPY))
    therapyGroup: dict[
        GkmTherapyGroupId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(TherapyGroup),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.THERAPY_GROUP))
    variantOrigin: dict[
        GkmVariantOriginId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(MappableConcept),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.VARIANT_ORIGIN))
    source: dict[
        GkmSourceId,
        Annotated[GkmBundleObject, _external_gkm_schema(Document)],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.SOURCE))
    method: dict[
        GkmMethodId,
        Annotated[GkmBundleObject, _external_gkm_schema(Method)],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.METHOD))
    organization: dict[
        GkmOrganizationId,
        Annotated[GkmBundleObject, _external_gkm_schema(Agent)],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.ORGANIZATION))
    proposition: dict[
        GkmPropositionId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(
                VariantClinicalSignificanceProposition,
                VariantDiagnosticProposition,
                VariantOncogenicityProposition,
                VariantPrognosticProposition,
                VariantTherapeuticResponseProposition,
            ),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.PROPOSITION))
    evidence: dict[
        GkmEvidenceId,
        Annotated[GkmBundleObject, _external_gkm_schema(Statement)],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.EVIDENCE))
    assertion: dict[
        GkmAssertionId,
        Annotated[
            GkmBundleObject,
            _external_gkm_schema(
                VariantClinicalSignificanceStatement,
                VariantOncogenicityStatement,
            ),
        ],
    ] = Field(json_schema_extra=_closed_collection_schema(Collection.ASSERTION))
    metadata: Metadata = Field(
        description="Bundle format, version, creation, and collection summary information."
    )
    failed_assertion_ids: list[int] = Field(
        description="CIViC Assertion IDs skipped because they could not be converted.",
    )
