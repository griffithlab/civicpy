"""Shared constants and enums for CIViC GKS exports."""

from enum import Enum


GA4GH_CURIE_PREFIX = "ga4gh"
_CIVIC_CURIE_PREFIX_ROOT = "civic"

ID_FIELD = "id"
CODE_FIELD = "code"
CODING_FIELD = "coding"
TYPE_FIELD = "type"
CONCEPT_TYPE_FIELD = "conceptType"
REFGET_ACCESSION_FIELD = "refgetAccession"
EXTENSIONS_FIELD = "extensions"
NAME_FIELD = "name"
DESCRIPTION_FIELD = "description"
ALIASES_FIELD = "aliases"
MAPPINGS_FIELD = "mappings"
CONDITIONS_FIELD = "conditions"
THERAPIES_FIELD = "therapies"
MEMBERSHIP_OPERATOR_FIELD = "membershipOperator"
PROPOSITION_FIELD = "proposition"
TARGET_PROPOSITION_FIELD = "targetProposition"
ALLELE_ORIGIN_QUALIFIER_FIELD = "alleleOriginQualifier"


class CivicGksBundleFormat(str, Enum):
    """Name and version of the CIViC GKS bundle format."""

    NAME = "civic-gks-bundle"
    VERSION = "0.1.0"


class CivicGksCuriePrefix(str, Enum):
    """CURIE prefixes used for CIViC objects represented in GKS.

    Each value is the prefix before the colon in a CIViC CURIE.
    """

    GENE = f"{_CIVIC_CURIE_PREFIX_ROOT}.gid"
    VARIANT = f"{_CIVIC_CURIE_PREFIX_ROOT}.vid"
    MOLECULAR_PROFILE = f"{_CIVIC_CURIE_PREFIX_ROOT}.mpid"
    DISEASE = f"{_CIVIC_CURIE_PREFIX_ROOT}.did"
    PHENOTYPE = f"{_CIVIC_CURIE_PREFIX_ROOT}.phenotype"
    THERAPY = f"{_CIVIC_CURIE_PREFIX_ROOT}.tid"
    SOURCE = f"{_CIVIC_CURIE_PREFIX_ROOT}.sid"
    EVIDENCE = f"{_CIVIC_CURIE_PREFIX_ROOT}.eid"
    ASSERTION = f"{_CIVIC_CURIE_PREFIX_ROOT}.aid"
    METHOD = f"{_CIVIC_CURIE_PREFIX_ROOT}.method"
    ORGANIZATION = f"{_CIVIC_CURIE_PREFIX_ROOT}.organization"
    VARIANT_ORIGIN = f"{_CIVIC_CURIE_PREFIX_ROOT}.variantOrigin"
    PROPOSITION = f"{_CIVIC_CURIE_PREFIX_ROOT}.proposition"
    CONDITION_SET = f"{_CIVIC_CURIE_PREFIX_ROOT}.conditionSet"
    THERAPY_GROUP = f"{_CIVIC_CURIE_PREFIX_ROOT}.therapyGroup"
