"""Module for representing CIViC assertion record as GKS AAC 2017 Study Statement"""

from abc import ABC
from enum import Enum
import re
from types import MappingProxyType

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import (
    Coding,
    ConceptMapping,
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
from pydantic import PrivateAttr
from civicpy.civic import (
    Assertion,
    Coordinate,
    Endorsement,
    Evidence,
    Disease,
    Organization,
    Therapy,
    GeneVariant,
    MolecularProfile,
    get_gene_by_id,
)


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

# SNP pattern
_SNP_RE = re.compile(r"RS\d+")


class _CivicGksAssertionRecord(ABC):
    """Abstract class for CIViC assertion record represented as GKS

    :param assertion: CIViC assertion record
    :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
        GKS object
    """

    _method: Method = PrivateAttr(
        Method(
            id="civic.method:2019",
            name="CIViC Curation SOP (2019)",
            reportedIn=Document(
                name="Danos et al., 2019, Genome Med.",
                title="Standard operating procedure for curation and clinical interpretation of variants in cancer",
                doi="10.1186/s13073-019-0687-x",
                pmid=31779674,
            ),
            methodType="variant curation standard operating procedure",
        )
    )
    _assertion: Assertion = PrivateAttr()

    def __init__(
        self,
        assertion: Assertion,
        endorsement: Endorsement | None = None,
    ) -> None:
        """Initialize _CivicGksAssertionRecord class

        :param assertion: CIViC assertion record
        :param endorsement: CIViC endorsement for the assertion, defaults to None
        :raises CivicGksRecordError: If CIViC assertion is not able to be represented as
            GKS object
        """
        object.__setattr__(self, "_assertion", assertion)

        if not self._assertion.is_valid_for_gks_json(emit_warnings=True):
            err_msg = "Assertion is not valid for GKS."
            raise CivicGksRecordError(err_msg)

        classification, strength = self.classification_and_strength(
            self._assertion.amp_level
        )

        if endorsement:
            organization = endorsement.organization
            contributions = [
                Contribution(
                    activityType=f"{endorsement.type}.last_reviewed",
                    date=endorsement.last_reviewed.split("T", 1)[0],
                    contributor=Agent(
                        id=f"civic.{organization.type}:{organization.id}",
                        name=organization.name,
                        description=organization.description,
                    ),
                )
            ]
        else:
            contributions = None

        super().__init__(
            id=f"civic.aid:{self._assertion.id}",
            contributions=contributions,
            description=self._assertion.description,
            specifiedBy=self.method(),
            proposition=self.proposition(self._assertion),
            direction=self.direction(self._assertion.assertion_direction),
            classification=classification,
            strength=strength,
            hasEvidenceLines=self.evidence_lines(),
        )

    def classification_and_strength(
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
            evidence_strength = self.evidence_strength(CivicEvidenceLevel(level))
            evidence_strength.primaryCoding.name = evidence_strength.name
            mappings = [
                ConceptMapping(
                    coding=evidence_strength.primaryCoding,
                    relation=Relation.EXACT_MATCH,
                )
            ]

            strength = MappableConcept(
                primaryCoding=Coding(code=Strength(f"Level {level}"), system=system),
                mappings=mappings,
            )

        return classification, strength

    def evidence_strength(self, evidence_level: CivicEvidenceLevel) -> MappableConcept:
        """Get CIViC Evidence Item strength

        :param evidence_level: CIViC evidence level
        :return: Strength for CIViC evidence item
        """
        return MappableConcept(
            name=CIVIC_EVIDENCE_LEVEL_TO_NAME[evidence_level],
            primaryCoding=Coding(
                system="https://civic.readthedocs.io/en/latest/model/evidence/level.html",
                code=evidence_level.value,
            ),
        )

    def method(self) -> Method:
        """Get GKS method

        :return: GKS method
        """
        return self._method.default

    def proposition(
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
        variant: GeneVariant = record.variants[0]

        params = {
            "subjectVariant": self.variant(record.molecular_profile, variant),
            "geneContextQualifier": self.gene(variant.gene.id),
            "alleleOriginQualifier": self.allele_origin_qualifier(record),
            "predicate": self.predicate(record),
        }

        record_type = (
            record.assertion_type
            if isinstance(record, Assertion)
            else record.evidence_type
        )

        if record_type == CivicEvidenceAssertionType.PREDICTIVE:
            condition_key = "conditionQualifier"
            params["objectTherapeutic"] = self.therapeutic(
                record.therapies, self._assertion.therapy_interaction_type
            )
            proposition = VariantTherapeuticResponseProposition
        else:
            condition_key = "objectCondition"

            if record_type == CivicEvidenceAssertionType.PROGNOSTIC:
                proposition = VariantPrognosticProposition
            else:
                proposition = VariantDiagnosticProposition

        params[condition_key] = self.condition(record.disease)

        return proposition(**params)

    def direction(self, record_direction: str) -> Direction | None:
        """Get direction for CIViC assertion or evidence item

        :param record_direction: CIViC assertion or evidence item's direction
        :return: Direction for CIViC assertion or evidence item
        """
        if record_direction == "SUPPORTS":
            return Direction.SUPPORTS
        if record_direction == "DOES_NOT_SUPPORT":
            return Direction.DISPUTES
        return None

    def evidence_lines(self) -> list[EvidenceLine]:
        """Get evidence lines for a CIViC assertion

        Only the CIViC evidence items that are supported for GKS will be included

        :raises CivicGksRecordError: If evidence item is not valid for GKS
        :return: List of CIViC evidence lines
        """
        evidence_lines: list[EvidenceLine] = []

        for evidence_item in self._assertion.evidence_items:
            if not evidence_item.is_valid_for_gks_json(emit_warnings=True):
                err_msg = f"Evidence {evidence_item.id} is not valid for GKS."
                raise CivicGksRecordError(err_msg)

            ev_source = evidence_item.source

            evidence_lines.append(
                EvidenceLine(
                    hasEvidenceItems=[
                        Statement(
                            id=f"civic.eid:{evidence_item.id}",
                            description=evidence_item.description,
                            specifiedBy=self.method(),
                            proposition=self.proposition(evidence_item),
                            direction=self.direction(evidence_item.evidence_direction),
                            strength=self.evidence_strength(
                                CivicEvidenceLevel(evidence_item.evidence_level)
                            ),
                            reportedIn=[
                                Document(
                                    id=f"civic.source:{ev_source.id}",
                                    name=ev_source.citation,
                                    title=ev_source.title,
                                    pmid=int(ev_source.citation_id)
                                    if ev_source.source_type == "PUBMED"
                                    else None,
                                )
                            ],
                        )
                    ],
                    directionOfEvidenceProvided=Direction.SUPPORTS,
                )
            )

        return evidence_lines

    @staticmethod
    def variant(
        molecular_profile: MolecularProfile, variant: GeneVariant
    ) -> CategoricalVariant:
        """Get GKS CategoricalVariant

        :param molecular_profile: CIViC simple molecular profile
        :param variant: CIViC variant corresponding to ``molecular_profile``
        :return: GKS CategoricalVariant
        """
        extensions = [
            Extension(
                name="CIViC Molecular Profile Score",
                value=molecular_profile.molecular_profile_score,
            )
        ]
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

        aliases = []
        mappings = []
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

        return CategoricalVariant(
            id=f"civic.mpid:{molecular_profile.id}",
            name=molecular_profile.name,
            description=molecular_profile.description,
            aliases=aliases or None,
            extensions=extensions,
            mappings=mappings or None,
        )

    @staticmethod
    def gene(gene_id: int) -> MappableConcept:
        """Get GKS gene

        :param gene_id: ID for CIViC gene
        :return: GKS gene
        """
        gene = get_gene_by_id(gene_id)

        if gene.aliases:
            extensions = [Extension(name="aliases", value=gene.aliases)]
        else:
            extensions = []

        if gene.description:
            extensions.append(Extension(name="description", value=gene.description))

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

        return MappableConcept(
            id=f"civic.gid:{gene_id}",
            conceptType="Gene",
            name=gene.name,
            mappings=mappings,
            extensions=extensions or None,
        )

    @staticmethod
    def allele_origin_qualifier(record: Evidence | Assertion) -> MappableConcept:
        """Get GKS allele origin qualifier

        :param record: CIViC assertion or evidence item
        :return: Allele origin qualifier
        """
        return MappableConcept(name=record.variant_origin)

    @staticmethod
    def predicate(
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
    def therapeutic(
        therapies: list[Therapy], therapy_interaction_type: str | None
    ) -> MappableConcept | TherapyGroup | None:
        """Get GKS therapeutic or TherapyGroup

        :param therapies: List of CIViC therapy records
        :param therapy_interaction_type: The interaction type for ``therapies``
        :raises CivicGksRecordError: If no therapies are found
        :return: GKS therapeutic (single therapy record) or TherapyGroup (multiple
            therapies)
        """

        def _get_therapy(therapy: Therapy) -> MappableConcept:
            """Get GKS therapeutic

            :param therapy: CIViC therapy
            :return: GKS therapeutic
            """
            if therapy.aliases:
                extensions = [Extension(name="aliases", value=therapy.aliases)]
            else:
                extensions = None

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

            return MappableConcept(
                id=f"civic.tid:{therapy.id}",
                name=therapy.name,
                conceptType="Therapy",
                mappings=mappings,
                extensions=extensions,
            )

        if not therapies:
            err_msg = "No therapies found"
            raise CivicGksRecordError(err_msg)

        if len(therapies) == 1:
            therapy: Therapy = therapies[0]
            return _get_therapy(therapy)
        else:
            membership_operator = (
                MembershipOperator.AND
                if therapy_interaction_type == CivicInteractionType.COMBINATION
                else MembershipOperator.OR
            )
            therapies_mc: list[MappableConcept] = [_get_therapy(t) for t in therapies]
            return TherapyGroup(
                therapies=therapies_mc, membershipOperator=membership_operator
            )

    @staticmethod
    def condition(disease: Disease) -> MappableConcept:
        """Get GKS condition

        :param disease: CIViC disease
        :return: GKS condition
        """
        if disease.doid:
            doid = f"DOID:{disease.doid}"
            mappings = [
                ConceptMapping(
                    coding=Coding(
                        id=doid,
                        code=doid,
                        system="https://disease-ontology.org/?id=",
                    ),
                    relation=Relation.EXACT_MATCH,
                )
            ]
        else:
            mappings = None

        return MappableConcept(
            id=f"civic.did:{disease.id}",
            conceptType="Disease",
            name=disease.name,
            mappings=mappings,
        )


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
