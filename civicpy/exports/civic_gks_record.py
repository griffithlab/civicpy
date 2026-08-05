"""Module for representing CIViC records as GKS representations

* CIViC Predictive, Prognostic, and Diagnostic Assertions map to Variant
  Clinical Significance Statements that follow the AMP/ASCO/CAP 2017 guidelines
* CIViC Oncogenic Assertions map to Variant Oncogenicity Statements that follow
  the ClinGen/CGC/VICC Oncogenicity 2022 guidelines
"""

import logging
import re
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType

from ga4gh.cat_vrs.models import (
    CategoricalVariant,
    Constraint,
    CopyChangeConstraint,
    DefiningAlleleConstraint,
    DefiningLocationConstraint,
    FeatureContextConstraint,
)
from ga4gh.cat_vrs.relations import LIFTOVER_TO_RELATION, TRANSLATION_OF_RELATION
from ga4gh.core.models import (
    Coding,
    ConceptMapping,
    Extension,
    MappableConcept,
    MembershipOperator,
    Relation,
    code,
    iriReference,
)
from ga4gh.va_spec.aac_2017 import (
    AMP_ASCO_CAP_CLASSIFICATION_MAP,
    AmpAscoCapClassificationCode,
    AmpAscoCapEvidenceLineStrength,
    DiagnosticEvidenceLine,
    PrognosticEvidenceLine,
    TherapeuticEvidenceLine,
    VariantClinicalSignificanceStatement,
)
from ga4gh.va_spec.base import (
    Agent,
    CcvClassification,
    ConditionSet,
    Contribution,
    DiagnosticPredicate,
    Direction,
    Document,
    Method,
    PrognosticPredicate,
    Statement,
    StrengthCode,
    System,
    TherapeuticResponsePredicate,
    TherapyGroup,
    VariantClinicalSignificanceProposition,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from ga4gh.va_spec.ccv_2022 import (
    METHOD as CCV_METHOD,
)
from ga4gh.va_spec.ccv_2022 import (
    VariantOncogenicityEvidenceLine,
    VariantOncogenicityStatement,
)
from ga4gh.va_spec.ccv_2022.derived_evidence import derive_onco_evidence_attributes
from ga4gh.vrs.models import Allele, CopyNumberChange, Expression, Syntax, Variation
from pydantic import BaseModel

from civicpy.civic import (
    LINKS_URL,
    Approval,
    Assertion,
    Coordinate,
    Disease,
    Evidence,
    Gene,
    GeneVariant,
    MolecularProfile,
    Organization,
    Phenotype,
    Source,
    Therapy,
)
from civicpy.exports.variation_normalizer import (
    VariationNormalizerDataProxy,
    VariationNormalizerRESTDataProxy,
)

_logger = logging.getLogger(__name__)

PUBMED_URL = "https://pubmed.ncbi.nlm.nih.gov"


class CivicGksRecordError(Exception):
    """Custom error for CIViC GKS Record exceptions"""


class CivicInteractionType(str, Enum):
    """Define constraints for the translation of supported CIViC interaction types into
    GKS

    SEQUENTIAL is not currently supported
    """

    SUBSTITUTES = "SUBSTITUTES"
    COMBINATION = "COMBINATION"


class CivicEvidenceAssertionType(str, Enum):
    """Define constraints for the translation of supported CIViC evidence and assertion
    types into GKS

    ONCOGENIC and PREDISPOSING are not currently supported
    """

    PREDICTIVE = "PREDICTIVE"
    PROGNOSTIC = "PROGNOSTIC"
    DIAGNOSTIC = "DIAGNOSTIC"
    ONCOGENIC = "ONCOGENIC"


class CivicSignificance(str, Enum):
    """Define constraints for significance values

    Not exhaustive. Only supports those that can be represented by GKS.
    """

    # Oncogenicity
    BENIGN = "BENIGN"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    LIKELY_ONCOGENIC = "LIKELY_ONCOGENIC"
    ONCOGENIC = "ONCOGENIC"
    UNCERTAIN_SIGNIFICANCE = "UNCERTAIN_SIGNIFICANCE"

    # Clinical Significance / impact
    BETTER_OUTCOME = "BETTER_OUTCOME"
    POOR_OUTCOME = "POOR_OUTCOME"
    POSITIVE = "POSITIVE"
    NEGATIVE = "NEGATIVE"
    RESISTANCE = "RESISTANCE"
    SENSITIVITY_RESPONSE = "SENSITIVITYRESPONSE"


CLINICAL_SIGNIFICANCE_ASSERTION_TYPES = [
    CivicEvidenceAssertionType.PREDICTIVE.value,
    CivicEvidenceAssertionType.PROGNOSTIC.value,
    CivicEvidenceAssertionType.DIAGNOSTIC.value,
]
ONCOGENIC_ASSERTION_TYPES = [CivicEvidenceAssertionType.ONCOGENIC.value]


class ClinVarSubmissionType(str, Enum):
    """Define supported submission types to ClinVar"""

    CLINICAL_IMPACT = "clinical_impact"
    ONCOGENICITY = "oncogenicity"


ASSERTION_TYPES_BY_CLINVAR_SUBMISSION_TYPE = MappingProxyType(
    {
        ClinVarSubmissionType.CLINICAL_IMPACT: CLINICAL_SIGNIFICANCE_ASSERTION_TYPES,
        ClinVarSubmissionType.ONCOGENICITY: ONCOGENIC_ASSERTION_TYPES,
    }
)


class CivicEvidenceLevel(str, Enum):
    """Define constraints for CIViC evidence levels"""

    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"


class CivicEvidenceName(str, Enum):
    """Define constraints for CIViC evidence names"""

    VALIDATED_ASSOCIATION = "Validated association"
    CLINICAL_EVIDENCE = "Clinical evidence"
    CASE_STUDY = "Case study"
    PRECLINICAL_EVIDENCE = "Preclinical evidence"
    INFERENTIAL_ASSOCIATION = "Inferential association"


# CIViC evidence level to CIViC evidence level name
CIVIC_EVIDENCE_LEVEL_TO_NAME = MappingProxyType(
    {
        CivicEvidenceLevel.A: CivicEvidenceName.VALIDATED_ASSOCIATION,
        CivicEvidenceLevel.B: CivicEvidenceName.CLINICAL_EVIDENCE,
        CivicEvidenceLevel.C: CivicEvidenceName.CASE_STUDY,
        CivicEvidenceLevel.D: CivicEvidenceName.PRECLINICAL_EVIDENCE,
        CivicEvidenceLevel.E: CivicEvidenceName.INFERENTIAL_ASSOCIATION,
    }
)

_IS_ONCOGENIC_FOR_PREDICATE = "isOncogenicFor"
# CIViC significance to GKS predicate
CLIN_SIG_TO_PREDICATE = MappingProxyType(
    {
        CivicSignificance.SENSITIVITY_RESPONSE.value: TherapeuticResponsePredicate.SENSITIVITY,
        CivicSignificance.RESISTANCE: TherapeuticResponsePredicate.RESISTANCE,
        CivicSignificance.POOR_OUTCOME: PrognosticPredicate.WORSE_OUTCOME,
        CivicSignificance.BETTER_OUTCOME: PrognosticPredicate.BETTER_OUTCOME,
        CivicSignificance.POSITIVE: DiagnosticPredicate.INCLUSIVE,
        CivicSignificance.NEGATIVE: DiagnosticPredicate.EXCLUSIVE,
        CivicSignificance.BENIGN: _IS_ONCOGENIC_FOR_PREDICATE,
        CivicSignificance.LIKELY_BENIGN: _IS_ONCOGENIC_FOR_PREDICATE,
        CivicSignificance.LIKELY_ONCOGENIC: _IS_ONCOGENIC_FOR_PREDICATE,
        CivicSignificance.ONCOGENIC: _IS_ONCOGENIC_FOR_PREDICATE,
        CivicSignificance.UNCERTAIN_SIGNIFICANCE: _IS_ONCOGENIC_FOR_PREDICATE,
    }
)

# CIViC variant origin to GKS allele origin (aligns with ClinVar API Schema)
VARIANT_ORIGIN_TO_ALLELE_ORIGIN = MappingProxyType(
    {
        "COMBINED": "unknown",
        "COMMON_GERMLINE": "germline",
        "MIXED": "unknown",
        "NA": "not applicable",
        "RARE_GERMLINE": "germline",
        "SOMATIC": "somatic",
        "UNKNOWN": "unknown",
    }
)


class CategoricalVariationType(str, Enum):
    """Supported categorical variation type extension values."""

    CATEGORICAL_CNV = "CategoricalCnv"
    FEATURE_CONTEXT = "FeatureContext"
    PROTEIN_SEQUENCE_CONSEQUENCE = "ProteinSequenceConsequence"
    UNDEFINED = "Undefined"


@dataclass(frozen=True)
class CategoricalVariantContext:
    """Container for Cat-VRS constraints and metadata."""

    constraints: list[Constraint]
    constraint_variation: Variation | None
    member_syntaxes: list[Syntax]
    categorical_variation_type: CategoricalVariationType


# SNP pattern
_SNP_RE = re.compile(r"RS\d+")


def resolve_variation_normalizer(
    variation_normalizer: VariationNormalizerDataProxy | None = None,
) -> VariationNormalizerDataProxy:
    """Return a provided VICC Variation Normalizer proxy or the default REST proxy.

    When generating multiple GKS records, callers should create one data proxy and
    pass it to each top-level GKS record. If omitted, each independently constructed
    top-level GKS record creates its own default REST proxy.

    :param variation_normalizer: VICC Variation Normalizer data proxy, if one is
        already configured for the export operation.
    :return: VICC Variation Normalizer data proxy to use.
    """
    if variation_normalizer is not None:
        return variation_normalizer

    return VariationNormalizerRESTDataProxy()


class CivicGksSop(Method):
    """Class for representing CIViC Curation SOP as GKS Method"""

    def __init__(self) -> None:
        """Initialize CivicGksSop class"""
        super().__init__(
            id="civic.method:2019",
            name="CIViC Curation SOP (2019)",
            reportedIn=Document(
                id="pmid:31779674",
                name="Danos et al., 2019, Genome Med.",
                title="Standard operating procedure for curation and clinical interpretation of variants in cancer",
                doi="10.1186/s13073-019-0687-x",
                pmid="31779674",
                urls=[
                    "https://doi.org/10.1186/s13073-019-0687-x",
                    f"{PUBMED_URL}/31779674/",
                ],
                aliases=["CIViC curation SOP"],
            ),
            methodType="curation",
        )


class CivicGksGene(MappableConcept):
    """Class for representing CIViC Gene as MappableConcept

    :param gene: CIViC gene record
    """

    def __init__(self, gene: Gene) -> None:
        """Initialize CivicGksGene class

        :param gene: CIViC gene record
        """
        super().__init__(
            id=f"civic.gid:{gene.id}",
            conceptType="Gene",
            name=gene.name,
            mappings=self.get_mappings(gene),
            extensions=self.get_extensions(gene),
        )

    @staticmethod
    def get_mappings(gene: Gene) -> list[ConceptMapping] | None:
        """Get mappings for CIViC gene

        :param gene: CIViC gene record
        :return: List of mappings containing entrez ID for CIViC gene, if found.
            Otherwise, ``None``.
        """
        if gene.entrez_id:
            entrez_id = str(gene.entrez_id)
            mappings = [
                ConceptMapping(
                    coding=Coding(
                        id=f"ncbigene:{entrez_id}",
                        code=entrez_id,
                        system="https://www.ncbi.nlm.nih.gov/gene/",
                    ),
                    relation=Relation.EXACT_MATCH,
                )
            ]
        else:
            mappings = None

        return mappings

    @staticmethod
    def get_extensions(gene: Gene) -> list[Extension] | None:
        """Get extensions for CIViC gene

        :param gene: CIViC gene record
        :return: List of extensions containing aliases and description for CIViC gene,
            if found. Otherwise, ``None``.
        """
        if gene.aliases:
            extensions = [Extension(name="aliases", value=gene.aliases)]
        else:
            extensions = []

        if gene.description:
            extensions.append(Extension(name="description", value=gene.description))

        return extensions or None


class CivicGksMolecularProfile(CategoricalVariant):
    """Represent a CIViC molecular profile as a Cat-VRS categorical variant.

    Simple molecular profiles containing one gene variant are represented as follows:

    * Gene-level mutation profiles use a feature-context constraint.
    * Profiles normalized to a VRS Allele use a defining-allele constraint and
      include successfully normalized coding and genomic HGVS expressions as
      members.
    * Profiles normalized to a VRS Copy Number Change use copy-change and
      defining-location constraints.
    * Profiles that cannot be normalized retain their CIViC metadata and use the
      ``Undefined`` categorical variation type without constraints or members.

    :param molecular_profile: CIViC molecular profile record
    :param variation_normalizer: VICC Variation Normalizer data proxy.
    """

    def __init__(
        self,
        molecular_profile: MolecularProfile,
        variation_normalizer: VariationNormalizerDataProxy | None = None,
    ) -> None:
        """Initialize CivicGksMolecularProfile class

        :param molecular_profile: CIViC molecular profile record
        :param variation_normalizer: VICC Variation Normalizer data proxy.
        :raises CivicGksRecordError: If molecular profile does not contain exactly one
            variant, or if the variant associated is not a Gene Variant
        """
        mp_id = molecular_profile.id
        mp_name = molecular_profile.sanitized_name
        variants = molecular_profile.variants

        if not variants:
            msg = f"Molecular profile contains no variants. mpid={mp_id}, name={mp_name!r}"
            raise CivicGksRecordError(msg)

        variant_count = len(variants)
        if variant_count != 1:
            msg = f"Only molecular profiles containing a single variant are supported. mpid={mp_id}, variant_count={variant_count}"
            raise CivicGksRecordError(msg)

        variant = molecular_profile.variants[0]
        if not isinstance(variant, GeneVariant):
            msg = f"Only GeneVariant records are supported. mpid={mp_id}, variant_type={type(variant).__name__!r}"
            raise CivicGksRecordError(msg)

        variation_normalizer = resolve_variation_normalizer(variation_normalizer)

        aliases, mappings = self._get_aliases_and_mappings(molecular_profile, variant)
        expressions = self._get_expressions(variant)
        extensions = self._get_extensions(molecular_profile, variant)

        categorical_variant_context = self._get_categorical_variant_context(
            molecular_profile, expressions, variation_normalizer
        )
        if categorical_variant_context:
            constraints = categorical_variant_context.constraints
            categorical_variation_type = (
                categorical_variant_context.categorical_variation_type
            )
            constraint_variation = categorical_variant_context.constraint_variation

            members = (
                self._build_members(
                    expressions,
                    categorical_variant_context.member_syntaxes,
                    constraint_variation,
                    variation_normalizer,
                )
                or None
            )
        else:
            constraints = None
            members = None
            categorical_variation_type = CategoricalVariationType.UNDEFINED

        extensions.append(
            Extension(
                name="categoricalVariationType",
                value=categorical_variation_type,
            )
        )

        super().__init__(
            id=f"civic.mpid:{molecular_profile.id}",
            name=molecular_profile.name,
            description=molecular_profile.description,
            aliases=aliases or None,
            extensions=extensions,
            mappings=mappings or None,
            constraints=constraints,
            members=members,
        )

    @staticmethod
    def _get_aliases_and_mappings(
        molecular_profile: MolecularProfile, variant: GeneVariant
    ) -> tuple[list[str], list[ConceptMapping]]:
        """Get aliases and mappings for a molecular profile

        :param molecular_profile: CIViC molecular profile record
        :param variant: Variant associated to molecular profile
        :return: A tuple containing aliases and dbSNP mappings for a molecular profile.
        """

        def get_variant_concept_mapping(variant: GeneVariant) -> ConceptMapping:
            """Build the CIViC concept mapping for a gene variant.

            :param variant: CIViC gene variant record.
            :return: Concept mapping for the gene variant.
            """
            extensions = [Extension(name="subtype", value=variant.subtype)]
            variant_types = [
                ConceptMapping(
                    coding=Coding(
                        id=f"civic.variant_type:{variant_type.id}",
                        code=variant_type.so_id,
                        name=variant_type.name,
                        system=f"{variant_type.url.rsplit('/', 1)[0]}/",
                    ),
                    relation=Relation.EXACT_MATCH,
                )
                for variant_type in variant.variant_types
                if variant_type.url is not None
            ]
            if variant_types:
                extensions.append(Extension(name="variant_types", value=variant_types))

            return ConceptMapping(
                coding=Coding(
                    id=f"civic.vid:{variant.id}",
                    code=str(variant.id),
                    name=variant.name,
                    system=f"{LINKS_URL}/variant/",
                    extensions=extensions,
                ),
                relation=Relation.EXACT_MATCH,
            )

        aliases = []
        mappings = [
            ConceptMapping(
                coding=Coding(
                    id=f"civic.mpid:{molecular_profile.id}",
                    code=str(molecular_profile.id),
                    system=f"{LINKS_URL}/molecular_profile/",
                ),
                relation=Relation.EXACT_MATCH,
            ),
            get_variant_concept_mapping(variant),
        ]

        allele_registry_id = variant.allele_registry_id
        if allele_registry_id:
            mappings.append(
                ConceptMapping(
                    coding=Coding(
                        id=f"clingen.allele:{allele_registry_id}",  # bioregistry
                        system="https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=",
                        code=allele_registry_id,
                    ),
                    relation=Relation.RELATED_MATCH,
                )
            )

        clinvar_ids = variant.clinvar_entries
        if clinvar_ids:
            mappings.extend(
                ConceptMapping(
                    coding=Coding(
                        id=f"clinvar:{clinvar_id}",  # identifiers.org + bioregistry
                        system="https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                        code=clinvar_id,
                    ),
                    relation=Relation.RELATED_MATCH,
                )
                for clinvar_id in clinvar_ids
                if clinvar_id and clinvar_id != "N/A"
            )

        for a in molecular_profile.aliases:
            if _SNP_RE.match(a):
                dbsnp_id = a.lower()
                mappings.append(
                    ConceptMapping(
                        coding=Coding(
                            id=f"dbsnp:{dbsnp_id}",  # identifiers.org + bioregistry
                            code=dbsnp_id,
                            system="https://www.ncbi.nlm.nih.gov/snp/",
                        ),
                        relation=Relation.RELATED_MATCH,
                    )
                )
            else:
                aliases.append(a)

        return aliases, mappings

    @staticmethod
    def _get_expressions(variant: GeneVariant) -> list[Expression]:
        """Get expressions for a variant

        :param variant: Variant associated to molecular profile
        :return: List of GKS expressions
        """

        def get_syntax(expr: str | None) -> Syntax | None:
            """Get syntax for an expression

            :param expr: HGVS expression
            :return: Syntax for HGVS expression, if p/c/g expression. Otherwise, None
            """
            if not expr:
                return

            if expr == "N/A":
                return

            if "p." in expr:
                return Syntax.HGVS_P

            if "c." in expr:
                return Syntax.HGVS_C

            if "g." in expr:
                return Syntax.HGVS_G

            return

        expressions_by_value: dict[str, Expression] = {}
        for expr, is_mane in [
            (item, False) for item in (variant.hgvs_expressions or [])
        ] + [(variant.mane_select_transcript, True)]:
            syntax = get_syntax(expr)
            if not syntax:
                continue

            extensions = (
                [Extension(name="isManeSelect", value=True)] if is_mane else None
            )
            expressions_by_value[expr] = Expression(
                syntax=syntax, value=expr, extensions=extensions
            )

        return list(expressions_by_value.values())

    @staticmethod
    def _get_extensions(
        molecular_profile: MolecularProfile,
        variant: GeneVariant,
    ) -> list[Extension]:
        """Get extensions for CIViC molecular profile

        :param molecular_profile: CIViC molecular profile record
        :param variant: Variant associated to molecular profile
        :return: List of extensions containing molecular profile score,
            HGVS descriptions, MANE Select transcript, and representative coordinate for a
            CIViC molecular profile record.
        """
        extensions = []

        for ext_name, ext_value in [
            (
                "CIViC Molecular Profile Score",
                molecular_profile.molecular_profile_score,
            ),
            ("hgvsDescriptions", variant.hgvs_expressions),
            ("maneSelectTranscript", variant.mane_select_transcript),
        ]:
            if ext_value is not None:
                extensions.append(Extension(name=ext_name, value=ext_value))

        if isinstance(variant.coordinates, Coordinate):
            coords = variant.coordinates
            ext_value = {
                "chromosome": coords.chromosome,
                "start": coords.start,
                "stop": coords.stop,
                "reference_bases": coords.reference_bases,
                "variant_bases": coords.variant_bases,
                "ensembl_version": coords.ensembl_version,
                "representative_transcript": coords.representative_transcript,
                "reference_build": coords.reference_build,
                "type": coords.type,
            }

            if all(v is not None for v in ext_value.values()):
                extensions.append(
                    Extension(
                        name="CIViC representative coordinate",
                        value=ext_value,
                    )
                )

        return extensions

    @staticmethod
    def _build_members(
        expressions: list[Expression],
        member_syntaxes: list[Syntax],
        constraint_variation: Variation | None,
        variation_normalizer: VariationNormalizerDataProxy,
    ) -> list[Variation]:
        """Build unique members from list of expressions.

        If multiple expressions normalize to the same VRS variation, they're
        merged into a single VRS variation. All expressions are retained and
        the lexicographically greatest expression is used as the name.

        The VRS variation used in a constraint is also included in the members.

        :param expressions: List of expressions for the gene variant
        :param member_syntaxes: Syntaxes that members will have
        :param constraint_variation: VRS variation used in a constraint.
        :param variation_normalizer: VICC Variation Normalizer data proxy.
        :return: Unique VRS members.
        """
        members_by_id = {}

        for expression in expressions:
            if expression.syntax not in member_syntaxes:
                continue

            hgvs_expr = expression.value
            normalized_variation = variation_normalizer.normalize(hgvs_expr)

            if not normalized_variation:
                continue

            vrs_variation = normalized_variation.model_copy(deep=True)
            variation_id = vrs_variation.id

            if variation_id in members_by_id:
                variation = members_by_id[variation_id]
                variation.root.expressions.append(expression)
                variation.root.name = max(variation.root.name, hgvs_expr)
            else:
                vrs_variation.name = hgvs_expr
                vrs_variation.expressions = [expression]
                members_by_id[variation_id] = Variation(root=vrs_variation)

        members = list(members_by_id.values())
        if constraint_variation:
            members.append(constraint_variation)

        return members

    def _get_categorical_variant_context(
        self,
        molecular_profile: MolecularProfile,
        expressions: list[Expression],
        variation_normalizer: VariationNormalizerDataProxy,
    ) -> CategoricalVariantContext | None:
        """Build the context needed for a categorical variant.

        :param molecular_profile: CIViC molecular profile record
        :param expressions: List of expressions for the gene variant
        :param variation_normalizer: VICC Variation Normalizer data proxy.
        :return: Categorical variant context, or ``None`` when the molecular profile
            cannot be normalized.
        :raises CivicGksRecordError: If a gene-level mutation has no gene mapping or
            the VICC Variation Normalizer returns an unsupported variation type.
        """
        parsed_name = molecular_profile.parsed_name
        if self._is_gene_level_mutation(parsed_name):
            return self._get_feature_context(parsed_name[0])

        normalized_variation = variation_normalizer.normalize_molecular_profile(
            molecular_profile
        )
        if not normalized_variation:
            _logger.warning(
                "Unable to normalize molecular profile to VRS variation. mpid=%s, name='%s'",
                molecular_profile.id,
                molecular_profile.name,
            )
            return None

        vrs_variation = normalized_variation.model_copy(deep=True)

        if isinstance(vrs_variation, Allele):
            return self._get_allele_context(vrs_variation, expressions)

        if isinstance(vrs_variation, CopyNumberChange):
            return self._get_copy_number_change_context(vrs_variation)

        msg = f"Unsupported VRS variation type returned by Variation Normalizer. mpid={molecular_profile.id}, type={type(vrs_variation).__name__!r}"
        raise CivicGksRecordError(msg)

    @staticmethod
    def _is_gene_level_mutation(parsed_name: list) -> bool:
        """Determine whether a parsed molecular profile is a gene-level mutation.

        :param parsed_name: Parsed molecular profile name components.
        :return: ``True`` if the profile represents a gene-level mutation.
        """
        return (
            len(parsed_name) == 2
            and isinstance(parsed_name[0], Gene)
            and parsed_name[1].name.lower() == "mutation"
        )

    @staticmethod
    def _get_feature_context(gene: Gene) -> CategoricalVariantContext:
        """Build Cat-VRS context for a gene-level mutation.

        :param gene: Gene associated with the molecular profile.
        :return: Categorical variant context containing a feature constraint.
        :raises CivicGksRecordError: If the gene has no concept mapping.
        """
        mappings = CivicGksGene.get_mappings(gene)
        if not mappings:
            msg = f"Unable to retrieve mappings for gene {gene.id}"
            raise CivicGksRecordError(msg)

        feature_context = FeatureContextConstraint(
            featureContext=MappableConcept(primaryCoding=mappings[0].coding)
        )
        return CategoricalVariantContext(
            constraints=[Constraint(root=feature_context)],
            member_syntaxes=[Syntax.HGVS_P, Syntax.HGVS_C, Syntax.HGVS_G],
            categorical_variation_type=CategoricalVariationType.FEATURE_CONTEXT,
            constraint_variation=None,
        )

    @staticmethod
    def _get_allele_context(
        allele: Allele, expressions: list[Expression]
    ) -> CategoricalVariantContext:
        """Build Cat-VRS context for a normalized allele.

        :param allele: Normalized VRS allele.
        :param expressions: HGVS expressions associated with the allele.
        :return: Categorical variant context containing a defining allele constraint.
        """
        protein_expressions = [
            expression
            for expression in expressions
            if expression.syntax == Syntax.HGVS_P
        ]
        allele.expressions = protein_expressions or None
        constraint = DefiningAlleleConstraint(
            allele=allele,
            relations=[LIFTOVER_TO_RELATION, TRANSLATION_OF_RELATION],
        )
        return CategoricalVariantContext(
            constraints=[Constraint(root=constraint)],
            member_syntaxes=[Syntax.HGVS_C, Syntax.HGVS_G],
            categorical_variation_type=CategoricalVariationType.PROTEIN_SEQUENCE_CONSEQUENCE,
            constraint_variation=Variation(root=allele),
        )

    @staticmethod
    def _get_copy_number_change_context(
        copy_number_change: CopyNumberChange,
    ) -> CategoricalVariantContext:
        """Build Cat-VRS context for a normalized copy-number change.

        :param copy_number_change: Normalized VRS copy-number change.
        :return: Categorical variant context containing copy and location constraints.
        """
        constraints = [
            Constraint(
                root=CopyChangeConstraint(copyChange=copy_number_change.copyChange)
            ),
            Constraint(
                root=DefiningLocationConstraint(
                    location=copy_number_change.location,
                    matchCharacteristic=MappableConcept(
                        primaryCoding=Coding(
                            code=code("is_within"),
                            system="ga4gh-gks-term:location-match",
                        )
                    ),
                    relations=[LIFTOVER_TO_RELATION],
                )
            ),
        ]
        return CategoricalVariantContext(
            constraints=constraints,
            member_syntaxes=[Syntax.HGVS_P, Syntax.HGVS_C, Syntax.HGVS_G],
            categorical_variation_type=CategoricalVariationType.CATEGORICAL_CNV,
            constraint_variation=Variation(root=copy_number_change),
        )


class CivicGksDisease(MappableConcept):
    """Class for representing CIViC Disease as MappableConcept

    :param disease: CIViC disease record
    """

    def __init__(self, disease: Disease) -> None:
        """Initialize CivicGksDisease class

        :param disease: CIViC disease record
        """
        super().__init__(
            id=f"civic.did:{disease.id}",
            conceptType="Disease",
            name=disease.name,
            mappings=self.get_mappings(disease),
        )

    @staticmethod
    def get_mappings(disease: Disease) -> list[ConceptMapping] | None:
        """Get mappings for CIViC disease

        :param disease: CIViC disease record
        :return: List of mappings containing DOID for CIViC disease, if found.
            Otherwise ``None``.
        """
        if disease.doid:
            mappings = [
                ConceptMapping(
                    coding=Coding(
                        code=f"DOID:{disease.doid}",
                        system="https://disease-ontology.org/?id=",
                    ),
                    relation=Relation.EXACT_MATCH,
                )
            ]
        else:
            mappings = None
        return mappings


class CivicGksPhenotype(MappableConcept):
    """Class for representing CIViC Phenotype as MappableConcept

    :param phenotype: CIViC phenotype record
    """

    def __init__(self, phenotype: Phenotype) -> None:
        """Initialize CivicGksPhenotype class

        :param phenotype: CIViC phenotype record
        """

        super().__init__(
            id=f"civic.{phenotype.type}:{phenotype.id}",
            conceptType=phenotype.type.capitalize(),
            name=phenotype.name,
            mappings=self.get_mappings(phenotype),
        )

    @staticmethod
    def get_mappings(phenotype: Phenotype) -> list[ConceptMapping]:
        """Get mappings for CIViC phenotype

        :param phenotype: phenotype disease record
        :return: List of mappings containing HPO ID for CIViC phenotype
        """
        _delimiter = "/"
        _system = phenotype.phenotype_url.rpartition(_delimiter)[0]
        return [
            ConceptMapping(
                coding=Coding(
                    code=phenotype.hpo_id,
                    system=f"{_system}{_delimiter}",
                ),
                relation=Relation.EXACT_MATCH,
            )
        ]


class CivicGksTherapy(MappableConcept):
    """Class for representing CIViC Therapy as MappableConcept

    :param therapy: CIViC therapy record
    """

    def __init__(self, therapy: Therapy) -> None:
        """Initialize CivicGksTherapy class

        :param therapy: CIViC therapy record
        """
        super().__init__(
            id=f"civic.tid:{therapy.id}",
            name=therapy.name,
            conceptType="Therapy",
            mappings=self.get_mappings(therapy),
            extensions=self.get_extensions(therapy),
        )

    @staticmethod
    def get_mappings(therapy: Therapy) -> list[ConceptMapping] | None:
        """Get mappings for CIViC therapy

        :param therapy: CIViC therapy record
        :return: List of mappings containing NCIt ID for CIViC therapy, if found.
            Otherwise ``None``.
        """
        if therapy.ncit_id:
            mappings = [
                ConceptMapping(
                    coding=Coding(
                        id=f"ncit:{therapy.ncit_id}",
                        code=therapy.ncit_id,
                        system="https://ncit.nci.nih.gov/ncitbrowser/ConceptReport.jsp?dictionary=NCI_Thesaurus&code=",
                    ),
                    relation=Relation.EXACT_MATCH,
                )
            ]
        else:
            mappings = None
        return mappings

    @staticmethod
    def get_extensions(therapy: Therapy) -> list[Extension] | None:
        """Get extensions for CIViC therapy

        :param therapy: CIViC therapy record
        :return: List of extensions containing aliases for a therapy.
        """
        if therapy.aliases:
            extensions = [Extension(name="aliases", value=therapy.aliases)]
        else:
            extensions = None

        return extensions


class CivicGksTherapyGroup(TherapyGroup):
    """Class for representing more than one CIViC therapies as a TherapyGroup

    :param therapies: List of CIViC therapy records
    :param therapy_interaction_type: Interaction type for list of therapies
    """

    def __init__(
        self, therapies: list[Therapy], therapy_interaction_type: str | None
    ) -> None:
        """Initialize CivicGksTherapyGroup class

        :param therapies: List of CIViC therapy records
        :param therapy_interaction_type: Interaction type for list of therapies
        :raises CivicGksRecordError: If no therapies were provided
        """
        if not therapies:
            err_msg = "No therapies provided"
            raise CivicGksRecordError(err_msg)

        membership_operator = (
            MembershipOperator.AND
            if therapy_interaction_type == CivicInteractionType.COMBINATION
            else MembershipOperator.OR
        )
        therapies_mc: list[MappableConcept] = [CivicGksTherapy(t) for t in therapies]

        super().__init__(therapies=therapies_mc, membershipOperator=membership_operator)


class _CivicGksEvidenceAssertionMixin:
    """Mixin for CIViC Evidence and Assertions"""

    @staticmethod
    def get_allele_origin_qualifier(record: Evidence | Assertion) -> MappableConcept:
        """Get GKS allele origin qualifier

        :param record: CIViC assertion or evidence item
        :return: Allele origin qualifier
        """
        variant_origin = record.variant_origin

        return MappableConcept(
            name=VARIANT_ORIGIN_TO_ALLELE_ORIGIN[variant_origin],
            mappings=[
                ConceptMapping(
                    coding=Coding(
                        code=variant_origin,
                        system="https://civicdb.org",
                        iris=[
                            iriReference(
                                root="https://civic.readthedocs.io/en/latest/model/evidence/origin.html"
                            )
                        ],
                    ),
                    relation=Relation.EXACT_MATCH,
                )
            ],
        )

    @staticmethod
    def get_predicate(
        record: Evidence | Assertion,
    ) -> (
        PrognosticPredicate
        | DiagnosticPredicate
        | TherapeuticResponsePredicate
        | str
        | None
    ):
        """Get GKS predicate

        :param record: CIViC assertion or evidence item
        :raises CivicGksRecordError: If significance is not supported for GKS
        :return: GKS predicate
        """
        try:
            return CLIN_SIG_TO_PREDICATE[record.significance]
        except KeyError:
            err_msg = f"Significance is not supported for GKS: {record.significance}"
            raise CivicGksRecordError(err_msg)

    @staticmethod
    def get_direction(record_direction: str) -> Direction | None:
        """Get direction for CIViC assertion or evidence item

        :param record_direction: CIViC assertion or evidence item's direction
        :return: Direction for CIViC assertion or evidence item
        """
        if record_direction == "SUPPORTS":
            return Direction.SUPPORTS
        if record_direction == "DOES_NOT_SUPPORT":
            return Direction.DISPUTES
        return None

    @staticmethod
    def get_evidence_strength(evidence_level: CivicEvidenceLevel) -> MappableConcept:
        """Get CIViC Evidence Item strength

        :param evidence_level: CIViC evidence level
        :return: Strength for CIViC evidence item
        """
        vicc_concept_vocab: ViccConceptVocab = VICC_CONCEPT_MAPPING[evidence_level]
        return MappableConcept(
            name=CIVIC_EVIDENCE_LEVEL_TO_NAME[evidence_level],
            primaryCoding=Coding(
                system="https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                code=evidence_level.value,
            ),
            mappings=[
                ConceptMapping(
                    coding=Coding(
                        system="https://go.osu.edu/evidence-codes",
                        code=vicc_concept_vocab.code,
                        name=vicc_concept_vocab.name,
                    ),
                    relation=Relation.EXACT_MATCH,
                )
            ],
        )

    def _get_proposition_params(
        self,
        record: Evidence | Assertion,
        record_type: CivicEvidenceAssertionType,
        variation_normalizer: VariationNormalizerDataProxy,
        is_clinical_significance_prop: bool = False,
    ) -> dict:
        """Get proposition parameters shared between propositions

        :param record: CIViC assertion or evidence item
        :param record_type: The type of ``record``
        :param variation_normalizer: Variation Normalizer data proxy
        :param is_clinical_significance_prop: Whether to generate parameters for
            VariantClinicalSignificanceProposition
        :return: Dictionary containing proposition parameters shared between
            propositions
        """
        variant: GeneVariant = record.molecular_profile.variants[0]

        params = {
            "subjectVariant": CivicGksMolecularProfile(
                record.molecular_profile, variation_normalizer
            ),
            "geneContextQualifier": CivicGksGene(variant.gene),
            "alleleOriginQualifier": self.get_allele_origin_qualifier(record),
            "predicate": self.get_predicate(record)
            if not is_clinical_significance_prop
            else VariantClinicalSignificanceProposition.model_fields[
                "predicate"
            ].default,
        }

        if record_type == CivicEvidenceAssertionType.ONCOGENIC:
            condition_key = "objectTumorType"
        elif (
            is_clinical_significance_prop
            or record_type != CivicEvidenceAssertionType.PREDICTIVE
        ):
            condition_key = "objectCondition"
        else:
            condition_key = "conditionQualifier"

        gks_disease = CivicGksDisease(record.disease)

        if record.phenotypes:
            conditions = [gks_disease]
            if len(record.phenotypes) > 1:
                conditions.append(
                    ConditionSet(
                        membershipOperator=MembershipOperator.OR,
                        conditions=[
                            CivicGksPhenotype(phenotype)
                            for phenotype in record.phenotypes
                        ],
                    )
                )
            else:
                conditions.append(CivicGksPhenotype(record.phenotypes[0]))

            params[condition_key] = ConditionSet(
                membershipOperator=MembershipOperator.AND, conditions=conditions
            )
        else:
            params[condition_key] = gks_disease
        return params

    def get_target_proposition(
        self,
        record: Evidence | Assertion,
        variation_normalizer: VariationNormalizerDataProxy,
    ) -> (
        VariantTherapeuticResponseProposition
        | VariantDiagnosticProposition
        | VariantPrognosticProposition
    ):
        """Get GKS target proposition

        :param record: CIViC assertion or evidence item
        :param variation_normalizer: Variation Normalizer data proxy
        :return: GKS target proposition
        """
        record_type = (
            record.assertion_type
            if isinstance(record, Assertion)
            else record.evidence_type
        )
        params: dict = self._get_proposition_params(
            record, record_type, variation_normalizer
        )

        if record_type == CivicEvidenceAssertionType.PREDICTIVE:
            if len(record.therapies) == 1:
                therapeutic = CivicGksTherapy(record.therapies[0])
            else:
                therapeutic = CivicGksTherapyGroup(
                    record.therapies, record.therapy_interaction_type
                )

            params["objectTherapeutic"] = therapeutic
            proposition_cls = VariantTherapeuticResponseProposition
        else:
            if record_type == CivicEvidenceAssertionType.PROGNOSTIC:
                proposition_cls = VariantPrognosticProposition
            else:
                proposition_cls = VariantDiagnosticProposition
        return proposition_cls(**params)


class CivicGksSource(Document):
    """Class for representing CIViC Source as Document

    :param source: CIViC source record
    """

    def __init__(self, source: Source, urls: list[str] | None = None) -> None:
        """Initialize CivicGksSource class

        :param source: CIViC source record
        :param urls: List of additional URLs to include in the document
        """
        source_urls = urls or []
        source_urls.extend([f"{LINKS_URL}/source/{source.id}", source.source_url])
        pmid = source.citation_id if source.source_type == "PUBMED" else None
        if pmc_id := source.pmc_id:
            source_urls.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}")

        super().__init__(
            id=f"civic.sid:{source.id}",
            name=source.citation,
            title=source.title,
            pmid=pmid,
            urls=source_urls,
        )


class ViccConceptVocab(BaseModel):
    """Define VICC Concept Vocab model (https://go.osu.edu/evidence-codes)"""

    code: str
    name: str


VICC_CONCEPT_MAPPING: MappingProxyType[CivicEvidenceLevel, ViccConceptVocab] = (
    MappingProxyType(
        {
            CivicEvidenceLevel.A: ViccConceptVocab(
                code="e000001", name="authoritative evidence"
            ),
            CivicEvidenceLevel.B: ViccConceptVocab(
                code="e000005", name="clinical cohort evidence"
            ),
            CivicEvidenceLevel.C: ViccConceptVocab(
                code="e000008", name="clinical case study evidence"
            ),
            CivicEvidenceLevel.D: ViccConceptVocab(
                code="e000009", name="preclinical evidence"
            ),
            CivicEvidenceLevel.E: ViccConceptVocab(
                code="e000010", name="inferential evidence"
            ),
        }
    )
)


class CivicGksEvidence(Statement, _CivicGksEvidenceAssertionMixin):
    """Class for representing CIViC Evidence item as Statement

    :param evidence_item: CIViC evidence item
    """

    def __init__(
        self,
        evidence_item: Evidence,
        variation_normalizer: VariationNormalizerDataProxy | None = None,
    ) -> None:
        """Initialize CivicGksEvidence class

        :param evidence_item: CIViC evidence item
        :param variation_normalizer: Variation Normalizer data proxy
        :raises CivicGksRecordError: If CIViC evidence item is not able to be
            represented as GKS object
        """
        if not evidence_item.is_valid_for_gks_json(emit_warnings=True):
            err_msg = f"Evidence {evidence_item.id} is not valid for GKS."
            raise CivicGksRecordError(err_msg)

        variation_normalizer = resolve_variation_normalizer(variation_normalizer)

        super().__init__(
            id=f"civic.eid:{evidence_item.id}",
            description=evidence_item.description,
            specifiedBy=CivicGksSop(),
            proposition=self.get_target_proposition(
                evidence_item, variation_normalizer
            ),
            direction=self.get_direction(evidence_item.evidence_direction),
            strength=self.get_evidence_strength(
                CivicEvidenceLevel(evidence_item.evidence_level)
            ),
            reportedIn=[
                CivicGksSource(
                    evidence_item.source,
                    urls=[f"{LINKS_URL}/evidence/{evidence_item.id}"],
                ),
            ],
        )


class _CivicGksAssertionMixin:
    """Mixin for CIViC Assertions"""

    @staticmethod
    def get_contributions(approval: Approval) -> list[Contribution]:
        """Get contributions for an approval

        :param approval: Approval for assertion
        :return: List of contributions, with one item containing when the approval was
            last reviewed an organization.
            Will include an extension, `is_approved_vcep`.
        """
        organization: Organization = approval.organization
        return [
            Contribution(
                activityType=f"{approval.type}.last_reviewed",
                date=approval.last_reviewed.split("T", 1)[0],
                contributor=Agent(
                    id=f"civic.{organization.type}:{organization.id}",
                    name=organization.name,
                    description=organization.description,
                    extensions=[
                        Extension(
                            name="is_approved_vcep", value=organization.is_approved_vcep
                        )
                    ],
                ),
            )
        ]

    @staticmethod
    def get_extensions(approval: Approval | None) -> list[Extension]:
        """Get extensions for an assertion

        :param approval: Approval for assertion, if exists
        :return: List of extensions for an assertion. This will contain a
            single record, `clinvar_accession` if one exists
        """
        extensions = []
        if approval and approval.clinvar_accession:
            extensions.append(
                Extension(name="clinvar_accession", value=approval.clinvar_accession)
            )
        return extensions

    @staticmethod
    def get_reported_in(assertion: Assertion) -> list[iriReference]:
        """Get reported in information for an assertion

        If multiple evidence items link to same source, will merge the source.

        :param assertion: CIViC assertion record
        :return: List of CIViC links to records which the assertion is reported in
        """
        reported_in: list[iriReference | Document] = [
            iriReference(f"{LINKS_URL}/assertion/{assertion.id}")
        ]
        civic_gks_sources = {}
        for evidence_item in assertion.evidence_items or []:
            source = evidence_item.source
            source_id = source.id
            evidence_item_url = f"{LINKS_URL}/evidence/{evidence_item.id}"

            if source_id in civic_gks_sources:
                civic_gks_sources[source_id].urls.append(evidence_item_url)
            else:
                civic_gks_sources[source_id] = Document.model_validate(
                    CivicGksSource(
                        source,
                        urls=[evidence_item_url],
                    )
                )
        reported_in.extend(list(civic_gks_sources.values()))

        return reported_in


class CivicGksClinSigAssertion(
    VariantClinicalSignificanceStatement,
    _CivicGksAssertionMixin,
    _CivicGksEvidenceAssertionMixin,
):
    """Class for CIViC predictive, prognostic, or diagnostic assertion record
    represented as GKS

    :param assertion: CIViC assertion record
    :param variation_normalizer: Variation Normalizer data proxy
    :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
        GKS object
    """

    def __init__(
        self,
        assertion: Assertion,
        approval: Approval | None = None,
        variation_normalizer: VariationNormalizerDataProxy | None = None,
    ) -> None:
        """Initialize CivicGksClinSigAssertion class

        :param assertion: CIViC assertion record
        :param approval: CIViC approval for the assertion, defaults to None
        :param variation_normalizer: Variation Normalizer data proxy
        :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
            GKS object
        """
        if assertion.assertion_type not in CLINICAL_SIGNIFICANCE_ASSERTION_TYPES:
            err_msg = (
                f"Assertion type must be one of {CLINICAL_SIGNIFICANCE_ASSERTION_TYPES}"
            )
            raise CivicGksRecordError(err_msg)

        if not assertion.is_valid_for_gks_json(emit_warnings=True):
            err_msg = "Assertion is not valid for GKS."
            raise CivicGksRecordError(err_msg)

        classification, strength, level = self.get_classification_strength_level(
            assertion.amp_level
        )
        contributions = self.get_contributions(approval) if approval else None
        variation_normalizer = resolve_variation_normalizer(variation_normalizer)

        super().__init__(
            id=f"civic.aid:{assertion.id}",
            contributions=contributions,
            description=assertion.description,
            specifiedBy=CivicGksSop(),
            proposition=self.get_proposition(assertion, variation_normalizer),
            direction=self.get_direction(assertion.assertion_direction),
            classification=classification,
            strength=strength,
            hasEvidenceLines=self.get_evidence_lines(
                assertion, level, variation_normalizer
            ),
            reportedIn=self.get_reported_in(assertion),
            extensions=self.get_extensions(approval) or None,
        )

    def get_classification_strength_level(
        self,
        amp_level: str,
    ) -> tuple[
        MappableConcept | None,
        MappableConcept | None,
        AmpAscoCapEvidenceLineStrength | None,
    ]:
        """Get classification, strength, and level

        :param amp_level: AMP/ASCO/CAP level
        :return: Classification, strength, and level, if found
        """
        classification = None
        strength = None
        system = System.AMP_ASCO_CAP
        level = None

        if amp_level != "NA":
            pattern = re.compile(r"TIER_(?P<tier>[IV]+)(?:_LEVEL_(?P<level>[A-D]))?")
            match = pattern.match(amp_level).groupdict()

            tier = AmpAscoCapClassificationCode(f"tier {match['tier'].lower()}")
            amp_asco_cap_config = AMP_ASCO_CAP_CLASSIFICATION_MAP[tier]

            classification = MappableConcept(
                name=amp_asco_cap_config.name,
                primaryCoding=Coding(code=tier, system=system),
            )

            strength = MappableConcept(
                primaryCoding=Coding(code=amp_asco_cap_config.strength, system=system)
            )
            level = AmpAscoCapEvidenceLineStrength(match["level"])

        return classification, strength, level

    def get_evidence_lines(
        self,
        assertion: Assertion,
        level: AmpAscoCapEvidenceLineStrength,
        variation_normalizer: VariationNormalizerDataProxy,
    ) -> (
        list[DiagnosticEvidenceLine]
        | list[PrognosticEvidenceLine]
        | list[TherapeuticEvidenceLine]
    ):
        """Get evidence lines for a CIViC assertion

        Only the CIViC evidence items that are supported for GKS will be included

        :param assertion: CIViC assertion
        :param level: The CIViC Assertion's AMP/ASCO/CAP category level
        :param variation_normalizer: Variation Normalizer data proxy
        :return: List of CIViC evidence lines
        :raises NotImplementedError: If the evidence line type is not supported.
        """
        direction = self.get_direction(assertion.assertion_direction)

        evidence_items: list[CivicGksEvidence] = []
        for evidence_item in assertion.evidence_items:
            try:
                evidence_items.append(
                    CivicGksEvidence(evidence_item, variation_normalizer)
                )
            except CivicGksRecordError as e:
                _logger.exception(
                    "Error translating %s to CivicGksEvidence: %s",
                    evidence_item.name,
                    str(e),
                )
            except Exception as e:
                _logger.exception(
                    "Unhandled error translating %s to CivicGksEvidence: %s",
                    evidence_item.name,
                    str(e),
                )

        if assertion.assertion_type == CivicEvidenceAssertionType.PREDICTIVE:
            evidence_line_cls = TherapeuticEvidenceLine
        elif assertion.assertion_type == CivicEvidenceAssertionType.DIAGNOSTIC:
            evidence_line_cls = DiagnosticEvidenceLine
        elif assertion.assertion_type == CivicEvidenceAssertionType.PROGNOSTIC:
            evidence_line_cls = PrognosticEvidenceLine
        else:
            msg = f"Evidence line type for assertion type is not supported: {assertion.assertion_type}"
            raise NotImplementedError(msg)

        return [
            evidence_line_cls(
                targetProposition=self.get_target_proposition(
                    assertion, variation_normalizer
                ),
                hasEvidenceItems=evidence_items or None,
                directionOfEvidenceProvided=direction,
                strengthOfEvidenceProvided=MappableConcept(
                    primaryCoding=(Coding(code=level, system=System.AMP_ASCO_CAP))
                ),
            ).root
        ]

    def get_proposition(
        self, assertion: Assertion, variation_normalizer: VariationNormalizerDataProxy
    ) -> VariantClinicalSignificanceProposition:
        """Get GKS proposition

        :param assertion: CIViC assertion record
        :param variation_normalizer: Variation Normalizer data proxy
        :return: GKS proposition
        """
        params = self._get_proposition_params(
            assertion,
            assertion.assertion_type,
            variation_normalizer,
            is_clinical_significance_prop=True,
        )
        return VariantClinicalSignificanceProposition(**params)


class CivicGksOncogenicAssertion(
    VariantOncogenicityStatement,
    _CivicGksAssertionMixin,
    _CivicGksEvidenceAssertionMixin,
):
    """Class for CIViC oncogenic assertion record represented as GKS"""

    def __init__(
        self,
        assertion: Assertion,
        approval: Approval | None = None,
        variation_normalizer: VariationNormalizerDataProxy | None = None,
    ) -> None:
        """Initialize CivicGksOncogenicAssertion class

        :param assertion: CIViC assertion record
        :param approval: CIViC approval for the assertion, defaults to None
        :param variation_normalizer: Variation Normalizer data proxy
        :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
            GKS object
        """
        if assertion.assertion_type not in ONCOGENIC_ASSERTION_TYPES:
            err_msg = f"Assertion type must be one of {ONCOGENIC_ASSERTION_TYPES}"
            raise CivicGksRecordError(err_msg)

        if not assertion.is_valid_for_gks_json(emit_warnings=True):
            err_msg = "Assertion is not valid for GKS."
            raise CivicGksRecordError(err_msg)

        variation_normalizer = resolve_variation_normalizer(variation_normalizer)

        contributions = self.get_contributions(approval) if approval else None
        proposition = self.get_proposition(assertion, variation_normalizer)
        classification, strength = self.get_classification_strength(
            assertion.significance
        )

        super().__init__(
            id=f"civic.aid:{assertion.id}",
            contributions=contributions,
            description=assertion.description,
            specifiedBy=CCV_METHOD,
            proposition=proposition,
            direction=self.get_direction(assertion.assertion_direction),
            classification=classification,
            strength=strength,
            hasEvidenceLines=self.get_evidence_lines(assertion),
            reportedIn=self.get_reported_in(assertion),
        )

    def get_classification_strength(
        self, significance
    ) -> tuple[MappableConcept, MappableConcept | None]:
        """Get classification and strength

        :param significance: Assertion's significance
        :return: Classification and strength, if found
        """
        _strength = None

        classification = MappableConcept(
            primaryCoding=Coding(
                code=code(CcvClassification[significance]), system=System.CCV
            )
        )

        if significance in {
            CivicSignificance.LIKELY_BENIGN,
            CivicSignificance.LIKELY_ONCOGENIC,
        }:
            _strength = StrengthCode.LIKELY
        elif significance in {CivicSignificance.BENIGN, CivicSignificance.ONCOGENIC}:
            _strength = StrengthCode.DEFINITIVE

        if _strength:
            strength = MappableConcept(
                primaryCoding=Coding(code=code(_strength.value), system=System.CCV)
            )
        else:
            strength = None

        return classification, strength

    def get_evidence_lines(
        self,
        assertion: Assertion,
    ) -> list[VariantOncogenicityEvidenceLine]:
        """Get evidence lines for a CIViC assertion

        :param assertion: CIViC assertion
        :return: List of CIViC evidence lines
        """
        direction = self.get_direction(assertion.assertion_direction)

        evidence_lines = []
        for clingen_code in assertion.clingen_codes or []:
            evidence_attrs = derive_onco_evidence_attributes(
                VariantOncogenicityEvidenceLine.Criterion(clingen_code.code)
            )
            evidence_lines.append(
                VariantOncogenicityEvidenceLine(
                    directionOfEvidenceProvided=direction,
                    **evidence_attrs.model_dump(),
                )
            )

        return evidence_lines

    def get_proposition(
        self, assertion: Assertion, variation_normalizer: VariationNormalizerDataProxy
    ) -> VariantOncogenicityProposition:
        """Get GKS proposition

        :param assertion: CIViC assertion record
        :param variation_normalizer: Variation Normalizer data proxy
        :return: GKS proposition
        """
        params = self._get_proposition_params(
            assertion,
            assertion.assertion_type,
            variation_normalizer,
            is_clinical_significance_prop=False,
        )
        return VariantOncogenicityProposition(**params)


def create_gks_record_from_assertion(
    assertion: Assertion,
    approval: Approval | None = None,
    submission_type_filter: ClinVarSubmissionType | None = None,
    variation_normalizer: VariationNormalizerDataProxy | None = None,
) -> CivicGksClinSigAssertion | CivicGksOncogenicAssertion:
    """Create GKS Record from CIViC Assertion

    :param assertion: CIViC assertion record
    :param approval: CIViC approval for the assertion, defaults to None
    :param submission_type_filter: Optional ClinVar submission type used to
        restrict which assertion types may be translated
    :param variation_normalizer: Variation Normalizer data proxy
    :raises NotImplementedError: If GKS Record translation is not yet supported.
        Currently, only the following assertion types are supported: DIAGNOSTIC,
        PREDICTIVE, PROGNOSTIC, and ONCOGENIC.
        Or if the assertion type is excluded by the provided ClinVar submission type
            filter.
    :return: GKS Assertion Record object
    """
    assertion_type = assertion.assertion_type

    if submission_type_filter:
        allowed_assertion_types = ASSERTION_TYPES_BY_CLINVAR_SUBMISSION_TYPE[
            submission_type_filter
        ]
        if assertion_type not in allowed_assertion_types:
            err_msg = f"Assertion type {assertion_type} is not supported for ClinVar submission type {submission_type_filter.value}"
            raise NotImplementedError(err_msg)

    if assertion_type in CLINICAL_SIGNIFICANCE_ASSERTION_TYPES:
        return CivicGksClinSigAssertion(
            assertion, approval=approval, variation_normalizer=variation_normalizer
        )

    if assertion_type in ONCOGENIC_ASSERTION_TYPES:
        return CivicGksOncogenicAssertion(
            assertion, approval=approval, variation_normalizer=variation_normalizer
        )

    err_msg = f"Assertion type {assertion_type} is not currently supported"
    raise NotImplementedError(err_msg)
