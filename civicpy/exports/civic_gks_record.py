"""Module for representing CIViC assertion record as GKS AAC 2017 Study Statement"""

from abc import ABC
from enum import Enum
import re
from types import MappingProxyType

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import (
    Coding,
    ConceptMapping,
    iriReference,
    Extension,
    MappableConcept,
    Relation,
)
from ga4gh.va_spec.aac_2017 import (
    Classification,
    Strength,
    VariantDiagnosticStudyStatement,
    VariantPrognosticStudyStatement,
    VariantTherapeuticResponseStudyStatement,
)
from ga4gh.va_spec.base import (
    Agent,
    Contribution,
    ConditionSet,
    DiagnosticPredicate,
    Direction,
    Document,
    EvidenceLine,
    MembershipOperator,
    Method,
    PrognosticPredicate,
    Statement,
    System,
    TherapeuticResponsePredicate,
    TherapyGroup,
    VariantDiagnosticProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from ga4gh.vrs.models import Expression, Syntax
from pydantic import BaseModel
from civicpy.civic import (
    LINKS_URL,
    Assertion,
    Coordinate,
    Approval,
    Evidence,
    Disease,
    Gene,
    Organization,
    Phenotype,
    Source,
    Therapy,
    GeneVariant,
    MolecularProfile,
)


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


# CIViC significance to GKS predicate
CLIN_SIG_TO_PREDICATE = MappingProxyType(
    {
        "SENSITIVITYRESPONSE": TherapeuticResponsePredicate.SENSITIVITY,
        "RESISTANCE": TherapeuticResponsePredicate.RESISTANCE,
        "POOR_OUTCOME": PrognosticPredicate.WORSE_OUTCOME,
        "BETTER_OUTCOME": PrognosticPredicate.BETTER_OUTCOME,
        "POSITIVE": DiagnosticPredicate.INCLUSIVE,
        "NEGATIVE": DiagnosticPredicate.EXCLUSIVE,
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


# SNP pattern
_SNP_RE = re.compile(r"RS\d+")


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

    def get_mappings(self, gene: Gene) -> list[ConceptMapping] | None:
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

    def get_extensions(self, gene: Gene) -> list[Extension] | None:
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
    """Class for representing CIViC Molecular Profile as CategoricalVariant

    :param molecular_profile: CIViC molecular profile record
    """

    def __init__(self, molecular_profile: MolecularProfile) -> None:
        """Initialize CivicGksMolecularProfile class

        :param molecular_profile: CIViC molecular profile record
        """
        aliases, mappings = self.get_aliases_and_mappings(molecular_profile)

        super().__init__(
            id=f"civic.mpid:{molecular_profile.id}",
            name=molecular_profile.name,
            description=molecular_profile.description,
            aliases=aliases or None,
            extensions=self.get_extensions(molecular_profile),
            mappings=mappings or None,
        )

    @staticmethod
    def get_aliases_and_mappings(
        molecular_profile: MolecularProfile,
    ) -> tuple[list[str], list[ConceptMapping]]:
        """Get aliases and mappings for a molecular profile

        :param molecular_profile: CIViC molecular profile record
        :return: A tuple containing aliases and dbSNP mappings for a molecular profile.
        """

        def _get_variant_concept_mapping(variant: GeneVariant) -> ConceptMapping:
            """Get concept mapping for variant

            :param variant: CIViC variant record
            :return: Concept mapping for CIViC variant record, containing variant
                subtype and variant types.
            """
            extensions = [Extension(name="subtype", value=variant.subtype)]

            variant_types = [
                ConceptMapping(
                    coding=Coding(
                        id=f"civic.variant_type:{vt.id}",
                        code=vt.so_id,
                        name=vt.name,
                        system=f"{vt.url.rsplit('/', 1)[0]}/",
                    ),
                    relation=Relation.EXACT_MATCH,
                )
                for vt in variant.variant_types
                if vt.url is not None
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
        variant: GeneVariant = molecular_profile.variants[0]
        variant_concept_mapping = _get_variant_concept_mapping(variant)
        mappings = [
            ConceptMapping(
                coding=Coding(
                    id=f"civic.mpid:{molecular_profile.id}",
                    code=str(molecular_profile.id),
                    system=f"{LINKS_URL}/molecular_profile/",
                ),
                relation=Relation.EXACT_MATCH,
            ),
            variant_concept_mapping,
        ]

        allele_registry_id = variant.allele_registry_id
        if allele_registry_id:
            mappings.append(
                ConceptMapping(
                    coding=Coding(
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
                        system="https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                        code=clinvar_id,
                    ),
                    relation=Relation.RELATED_MATCH,
                )
                for clinvar_id in clinvar_ids
            )

        for a in molecular_profile.aliases:
            if _SNP_RE.match(a):
                a = a.lower()
                mappings.append(
                    ConceptMapping(
                        coding=Coding(
                            code=a,
                            system="https://www.ncbi.nlm.nih.gov/snp/",
                        ),
                        relation=Relation.RELATED_MATCH,
                    )
                )
            else:
                aliases.append(a)

        return aliases, mappings

    @staticmethod
    def get_extensions(molecular_profile: MolecularProfile) -> list[Extension]:
        """Get extensions for CIViC molecular profile

        :param molecular_profile: CIViC molecular profile record
        :return: List of extensions containing molecular profile score, expressions,
            and representative for a CIViC molecular profile record.
        """
        extensions = [
            Extension(
                name="CIViC Molecular Profile Score",
                value=molecular_profile.molecular_profile_score,
            )
        ]

        variant: GeneVariant = molecular_profile.variants[0]
        if variant.hgvs_expressions:
            expressions = []

            for hgvs_expr in variant.hgvs_expressions:
                if hgvs_expr == "N/A":
                    continue

                if "p." in hgvs_expr:
                    syntax = Syntax.HGVS_P
                elif "c." in hgvs_expr:
                    syntax = Syntax.HGVS_C
                elif "g." in hgvs_expr:
                    syntax = Syntax.HGVS_G
                else:
                    continue

                expressions.append(Expression(syntax=syntax, value=hgvs_expr))

            if expressions:
                extensions.append(Extension(name="expressions", value=expressions))

        if isinstance(variant.coordinates, Coordinate):
            coords = variant.coordinates
            extensions.append(
                Extension(
                    name="CIViC representative coordinate",
                    value={
                        "chromosome": coords.chromosome,
                        "start": coords.start,
                        "stop": coords.stop,
                        "reference_bases": coords.reference_bases,
                        "variant_bases": coords.variant_bases,
                        "ensembl_version": coords.ensembl_version,
                        "representative_transcript": coords.representative_transcript,
                        "reference_build": coords.reference_build,
                        "type": coords.type,
                    },
                )
            )
        return extensions


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
    @staticmethod
    def get_allele_origin_qualifier(record: Evidence | Assertion) -> MappableConcept:
        """Get GKS allele origin qualifier

        :param record: CIViC assertion or evidence item
        :return: Allele origin qualifier
        """
        variant_origin = record.variant_origin

        return MappableConcept(
            name=VARIANT_ORIGIN_TO_ALLELE_ORIGIN[variant_origin],
            extensions=[Extension(name="civic_variant_origin", value=variant_origin)],
        )

    @staticmethod
    def get_predicate(
        record: Evidence | Assertion,
    ) -> (
        PrognosticPredicate | DiagnosticPredicate | TherapeuticResponsePredicate | None
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

    def get_proposition(
        self, record: Evidence | Assertion
    ) -> (
        VariantTherapeuticResponseProposition
        | VariantDiagnosticProposition
        | VariantPrognosticProposition
    ):
        """Get GKS proposition

        :param record: CIViC assertion or evidence item
        :return: GKS proposition
        """
        variant: GeneVariant = record.molecular_profile.variants[0]

        params = {
            "subjectVariant": CivicGksMolecularProfile(record.molecular_profile),
            "geneContextQualifier": CivicGksGene(variant.gene),
            "alleleOriginQualifier": self.get_allele_origin_qualifier(record),
            "predicate": self.get_predicate(record),
        }

        record_type = (
            record.assertion_type
            if isinstance(record, Assertion)
            else record.evidence_type
        )

        if record_type == CivicEvidenceAssertionType.PREDICTIVE:
            condition_key = "conditionQualifier"
            if len(record.therapies) == 1:
                therapeutic = CivicGksTherapy(record.therapies[0])
            else:
                therapeutic = CivicGksTherapyGroup(
                    record.therapies, record.therapy_interaction_type
                )

            params["objectTherapeutic"] = therapeutic
            proposition = VariantTherapeuticResponseProposition
        else:
            condition_key = "objectCondition"

            if record_type == CivicEvidenceAssertionType.PROGNOSTIC:
                proposition = VariantPrognosticProposition
            else:
                proposition = VariantDiagnosticProposition

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
        return proposition(**params)


class CivicGksSource(Document):
    """Class for representing CIViC Source as Document

    :param source: CIViC source record
    """

    def __init__(self, source: Source) -> None:
        """Initialize CivicGksSource class

        :param source: CIViC source record
        """
        urls = [f"{LINKS_URL}/source/{source.id}", source.source_url]
        pmid = source.citation_id if source.source_type == "PUBMED" else None
        if pmc_id := source.pmc_id:
            urls.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}")

        super().__init__(
            id=f"civic.sid:{source.id}",
            name=source.citation,
            title=source.title,
            pmid=pmid,
            urls=urls,
        )


class ViccConceptVocab(BaseModel):
    """Define VICC Concept Vocab model (https://go.osu.edu/evidence-codes)"""

    code: str
    name: str


VICC_CONCEPT_MAPPING: dict[CivicEvidenceLevel, ViccConceptVocab] = MappingProxyType(
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


class CivicGksEvidence(Statement, _CivicGksEvidenceAssertionMixin):
    """Class for representing CIViC Evidence item as Statement

    :param evidence_item: CIViC evidence item
    """

    def __init__(self, evidence_item: Evidence) -> None:
        """Initialize CivicGksEvidence class

        :param evidence_item: CIViC evidence item
        :raises CivicGksRecordError: If CIViC evidence item is not able to be
            represented as GKS object
        """
        if not evidence_item.is_valid_for_gks_json(emit_warnings=True):
            err_msg = f"Evidence {evidence_item.id} is not valid for GKS."
            raise CivicGksRecordError(err_msg)

        super().__init__(
            id=f"civic.eid:{evidence_item.id}",
            description=evidence_item.description,
            specifiedBy=CivicGksSop(),
            proposition=self.get_proposition(evidence_item),
            direction=self.get_direction(evidence_item.evidence_direction),
            strength=self.get_evidence_strength(
                CivicEvidenceLevel(evidence_item.evidence_level)
            ),
            reportedIn=[
                CivicGksSource(evidence_item.source),
                iriReference(f"{LINKS_URL}/evidence/{evidence_item.id}"),
            ],
        )


class _CivicGksAssertionRecord(_CivicGksEvidenceAssertionMixin, ABC):
    """Abstract class for CIViC assertion record represented as GKS

    :param assertion: CIViC assertion record
    :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
        GKS object
    """

    def __init__(
        self,
        assertion: Assertion,
        approval: Approval | None = None,
    ) -> None:
        """Initialize _CivicGksAssertionRecord class

        :param assertion: CIViC assertion record
        :param approval: CIViC approval for the assertion, defaults to None
        :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
            GKS object
        """
        if not assertion.is_valid_for_gks_json(emit_warnings=True):
            err_msg = "Assertion is not valid for GKS."
            raise CivicGksRecordError(err_msg)

        classification, strength = self.get_classification_and_strength(
            assertion.amp_level
        )
        contributions = self.get_contributions(approval) if approval else None

        super().__init__(
            id=f"civic.aid:{assertion.id}",
            contributions=contributions,
            description=assertion.description,
            specifiedBy=CivicGksSop(),
            proposition=self.get_proposition(assertion),
            direction=self.get_direction(assertion.assertion_direction),
            classification=classification,
            strength=strength,
            hasEvidenceLines=self.get_evidence_lines(assertion, strength),
            reportedIn=[iriReference(f"{LINKS_URL}/assertion/{assertion.id}")],
        )

    @staticmethod
    def get_contributions(approval: Approval) -> list[Contribution]:
        """Get contributions for an approval

        :param approval: Approval for assertion
        :return: List of contributions containing when the approval was last reviewed
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
                    extensions=[Extension(name="is_approved_vcep", value=organization.is_approved_vcep)],
                ),
            )
        ]

    def get_classification_and_strength(
        self,
        amp_level: str,
    ) -> tuple[MappableConcept | None, MappableConcept | None]:
        """Get classification and strength

        :param amp_level: AMP/ASCO/CAP level
        :return: Classification and strength, if found
        """
        classification = None
        strength = None
        system = System.AMP_ASCO_CAP

        if amp_level != "NA":
            pattern = re.compile(r"TIER_(?P<tier>[IV]+)(?:_LEVEL_(?P<level>[A-D]))?")
            match = pattern.match(amp_level).groupdict()
            classification = MappableConcept(
                primaryCoding=Coding(
                    code=Classification(f"Tier {match['tier']}"), system=system
                ),
            )

            level = match["level"]
            strength = MappableConcept(
                primaryCoding=Coding(code=Strength(f"Level {level}"), system=system)
            )

        return classification, strength

    def get_evidence_lines(
        self, assertion: Assertion, strength: MappableConcept
    ) -> list[EvidenceLine]:
        """Get evidence lines for a CIViC assertion

        Only the CIViC evidence items that are supported for GKS will be included

        :param assertion: CIViC assertion
        :param strength: The CIViC Assertion's strength
        :return: List of CIViC evidence lines
        """
        direction = (
            Direction.SUPPORTS
            if assertion.assertion_direction == "SUPPORTS"
            else Direction.DISPUTES
        )
        evidence_items = [
            CivicGksEvidence(evidence_item)
            for evidence_item in assertion.evidence_items
        ]

        return [
            EvidenceLine(
                hasEvidenceItems=evidence_items,
                directionOfEvidenceProvided=direction,
                strengthOfEvidenceProvided=strength,
            )
        ]


class CivicGksPredictiveAssertion(
    _CivicGksAssertionRecord, VariantTherapeuticResponseStudyStatement
):
    """Class for representing CIViC predictive assertion as GKS VariantTherapeuticResponseStudyStatement"""


class CivicGksDiagnosticAssertion(
    _CivicGksAssertionRecord, VariantDiagnosticStudyStatement
):
    """Class for representing CIViC diagnostic assertion as GKS VariantDiagnosticStudyStatement"""


class CivicGksPrognosticAssertion(
    _CivicGksAssertionRecord, VariantPrognosticStudyStatement
):
    """Class for representing CIViC prognostic assertion as GKS VariantPrognosticStudyStatement"""


def create_gks_record_from_assertion(
    assertion: Assertion, approval: Approval | None = None
) -> (
    CivicGksDiagnosticAssertion
    | CivicGksPredictiveAssertion
    | CivicGksPrognosticAssertion
):
    """Create GKS Record from CIViC Assertion

    :param assertion: CIViC assertion record
    :param approval: CIViC approval for the assertion, defaults to None
    :raises NotImplementedError: If GKS Record translation is not yet supported.
        Currently, only the following assertion types are supported: DIAGNOSTIC,
        PREDICTIVE, and PROGNOSTIC.
    :return: GKS Assertion Record object
    """
    if assertion.assertion_type == "DIAGNOSTIC":
        return CivicGksDiagnosticAssertion(assertion, approval=approval)

    if assertion.assertion_type == "PREDICTIVE":
        return CivicGksPredictiveAssertion(assertion, approval=approval)

    if assertion.assertion_type == "PROGNOSTIC":
        return CivicGksPrognosticAssertion(assertion, approval=approval)

    err_msg = f"Assertion type {assertion.assertion_type} is not currently supported"
    raise NotImplementedError(err_msg)
