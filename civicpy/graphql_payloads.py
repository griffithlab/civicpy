def _construct_get_gene_payload():
    return """
        query gene($id: Int!) {
            gene(id: $id) {
                id
                name
                description
                entrez_id: entrezId
                aliases: featureAliases
                sources {
                    id
                }
            }
        }"""


def _construct_get_all_genes_payload():
    return """
        query genes($after: String) {
            genes(after: $after, evidenceStatusFilter: ALL) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                description
                entrez_id: entrezId
                aliases: featureAliases
                sources {
                    id
                }
              }
            }
        }"""


def _construct_get_factor_payload():
    return """
        query factor($id: Int!) {
            factor(id: $id) {
                name
                full_name: fullName
                description
                ncit_id: ncitId
                aliases: featureAliases
                sources {
                    id
                }
            }
        }"""


def _construct_get_all_factors_payload():
    return """
        query factors($after: String) {
            factors(after: $after, evidenceStatusFilter: ALL) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                full_name: fullName
                description
                ncit_id: ncitId
                aliases: featureAliases
                sources {
                    id
                }
              }
            }
        }"""


def _construct_get_fusion_payload():
    return """
        query fusion($id: Int!) {
            fusion(id: $id) {
                id
                name
                description
                threePrimeGene {
                    id
                }
                fivePrimeGene {
                    id
                }
                three_prime_partner_status: threePrimePartnerStatus
                five_prime_partner_status: fivePrimePartnerStatus
                aliases: featureAliases
                sources {
                    id
                }
            }
        }"""


def _construct_get_all_fusions_payload():
    return """
        query fusions($after: String) {
            fusions(after: $after, evidenceStatusFilter: ALL) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                description
                threePrimeGene {
                    id
                }
                fivePrimeGene {
                    id
                }
                three_prime_partner_status: threePrimePartnerStatus
                five_prime_partner_status: fivePrimePartnerStatus
                aliases: featureAliases
                sources {
                    id
                }
              }
            }
        }"""


def _construct_get_molecular_profile_payload():
    return """
        query molecularProfile($id: Int!) {
            molecular_profile: molecularProfile(id: $id) {
                id
                description
                molecular_profile_score: molecularProfileScore
                name
                variants {
                  id
                }
                aliases: molecularProfileAliases
                parsed_name: parsedName {
                    type: __typename
                    ... on MolecularProfileTextSegment {
                        text
                    }
                    ... on Feature {
                        id
                        name
                        featureType
                    }
                    ... on Variant {
                        id
                        name
                        deprecated
                    }
                }
                sources {
                    id
                }
            }
        }"""



def _construct_get_all_molecular_profiles_payload():
    return """
        query molecularProfiles($after: String) {
            molecular_profiles: molecularProfiles(after: $after, evidenceStatusFilter: ALL) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                description
                molecular_profile_score: molecularProfileScore
                name
                variants {
                  id
                }
                aliases: molecularProfileAliases
                parsed_name: parsedName {
                    type: __typename
                    ... on MolecularProfileTextSegment {
                        text
                    }
                    ... on Feature {
                        id
                        name
                        featureType
                    }
                    ... on Variant {
                        id
                        name
                        deprecated
                    }
                }
                sources {
                    id
                }
              }
            }
        }"""


def _construct_get_variant_payload():
    return """
        query variant($id: Int!) {
            variant(id: $id) {
                __typename
                id
                name
                ... on GeneVariant {
                    allele_registry_id: alleleRegistryId
                    clinvar_entries: clinvarIds
                    hgvs_expressions: hgvsDescriptions
                    coordinates {
                        reference_build: referenceBuild
                        ensembl_version: ensemblVersion
                        chromosome
                        representative_transcript: representativeTranscript
                        start
                        stop
                        reference_bases: referenceBases
                        variant_bases: variantBases
                    }
                }
                ... on FusionVariant {
                    vicc_compliant_name: viccCompliantName
                    five_prime_coordinates: fivePrimeCoordinates {
                        reference_build: referenceBuild
                        ensembl_version: ensemblVersion
                        chromosome
                        representative_transcript: representativeTranscript
                        start
                        stop
                        reference_bases: referenceBases
                        variant_bases: variantBases
                    }
                    three_prime_coordinates: threePrimeCoordinates {
                        reference_build: referenceBuild
                        ensembl_version: ensemblVersion
                        chromosome
                        representative_transcript: representativeTranscript
                        start
                        stop
                        reference_bases: referenceBases
                        variant_bases: variantBases
                    }
                    five_prime_start_exon_coordinates: fivePrimeStartExonCoordinates {
                        chromosome
                        ensembl_id: ensemblId
                        ensembl_version: ensemblVersion
                        exon
                        exon_offset: exonOffset
                        exon_offset_direction: exonOffsetDirection
                        reference_build: referenceBuild
                        representative_transcript: representativeTranscript
                        start
                        stop
                        strand
                    }
                    five_prime_end_exon_coordinates: fivePrimeEndExonCoordinates {
                        chromosome
                        ensembl_id: ensemblId
                        ensembl_version: ensemblVersion
                        exon
                        exon_offset: exonOffset
                        exon_offset_direction: exonOffsetDirection
                        reference_build: referenceBuild
                        representative_transcript: representativeTranscript
                        start
                        stop
                        strand
                    }
                    three_prime_start_exon_coordinates: threePrimeStartExonCoordinates {
                        chromosome
                        ensembl_id: ensemblId
                        ensembl_version: ensemblVersion
                        exon
                        exon_offset: exonOffset
                        exon_offset_direction: exonOffsetDirection
                        reference_build: referenceBuild
                        representative_transcript: representativeTranscript
                        start
                        stop
                        strand
                    }
                    three_prime_end_exon_coordinates: threePrimeEndExonCoordinates {
                        chromosome
                        ensembl_id: ensemblId
                        ensembl_version: ensemblVersion
                        exon
                        exon_offset: exonOffset
                        exon_offset_direction: exonOffsetDirection
                        reference_build: referenceBuild
                        representative_transcript: representativeTranscript
                        start
                        stop
                        strand
                    }
                }
                ... on FactorVariant {
                    ncit_id: ncitId
                }
                feature {
                    id
                    name
                    featureInstance {
                        ... on Gene {
                            entrezId
                        }
                    }
                }
                single_variant_molecular_profile_id: singleVariantMolecularProfileId
                variant_aliases: variantAliases
                variant_types: variantTypes {
                    id
                    name
                    so_id: soid
                    description
                    url
                }
            }
        }"""


def _construct_get_all_variants_payload():
    return """
        query variants($after: String) {
            variants(after: $after) {
                totalCount
                pageInfo {
                  hasNextPage
                  endCursor
                }
                nodes {
                    __typename
                    id
                    name
                    ... on GeneVariant {
                        allele_registry_id: alleleRegistryId
                        clinvar_entries: clinvarIds
                        hgvs_expressions: hgvsDescriptions
                        coordinates {
                            reference_build: referenceBuild
                            ensembl_version: ensemblVersion
                            chromosome
                            representative_transcript: representativeTranscript
                            start
                            stop
                            reference_bases: referenceBases
                            variant_bases: variantBases
                        }
                    }
                    ... on FactorVariant {
                        ncit_id: ncitId
                    }
                    ... on FusionVariant {
                        vicc_compliant_name: viccCompliantName
                        five_prime_coordinates: fivePrimeCoordinates {
                            reference_build: referenceBuild
                            ensembl_version: ensemblVersion
                            chromosome
                            representative_transcript: representativeTranscript
                            start
                            stop
                            reference_bases: referenceBases
                            variant_bases: variantBases
                        }
                        three_prime_coordinates: threePrimeCoordinates {
                            reference_build: referenceBuild
                            ensembl_version: ensemblVersion
                            chromosome
                            representative_transcript: representativeTranscript
                            start
                            stop
                            reference_bases: referenceBases
                            variant_bases: variantBases
                        }
                        five_prime_start_exon_coordinates: fivePrimeStartExonCoordinates {
                            chromosome
                            ensembl_id: ensemblId
                            ensembl_version: ensemblVersion
                            exon
                            exon_offset: exonOffset
                            exon_offset_direction: exonOffsetDirection
                            reference_build: referenceBuild
                            representative_transcript: representativeTranscript
                            start
                            stop
                            strand
                        }
                        five_prime_end_exon_coordinates: fivePrimeEndExonCoordinates {
                            chromosome
                            ensembl_id: ensemblId
                            ensembl_version: ensemblVersion
                            exon
                            exon_offset: exonOffset
                            exon_offset_direction: exonOffsetDirection
                            reference_build: referenceBuild
                            representative_transcript: representativeTranscript
                            start
                            stop
                            strand
                        }
                        three_prime_start_exon_coordinates: threePrimeStartExonCoordinates {
                            chromosome
                            ensembl_id: ensemblId
                            ensembl_version: ensemblVersion
                            exon
                            exon_offset: exonOffset
                            exon_offset_direction: exonOffsetDirection
                            reference_build: referenceBuild
                            representative_transcript: representativeTranscript
                            start
                            stop
                            strand
                        }
                        three_prime_end_exon_coordinates: threePrimeEndExonCoordinates {
                            chromosome
                            ensembl_id: ensemblId
                            ensembl_version: ensemblVersion
                            exon
                            exon_offset: exonOffset
                            exon_offset_direction: exonOffsetDirection
                            reference_build: referenceBuild
                            representative_transcript: representativeTranscript
                            start
                            stop
                            strand
                        }
                    }
                    feature {
                        id
                        name
                        featureInstance {
                            ... on Gene {
                                entrezId
                            }
                        }
                    }
                    single_variant_molecular_profile_id: singleVariantMolecularProfileId
                    variant_aliases: variantAliases
                    variant_types: variantTypes {
                        id
                        name
                        so_id: soid
                        description
                        url
                    }
                }
            }
        }"""

def _construct_get_evidence_payload():
    return """
        query evidenceItem($id: Int!) {
            evidence: evidenceItem(id: $id) {
                id
                name
                significance
                description
                therapy_interaction_type: therapyInteractionType
                evidence_direction: evidenceDirection
                evidence_level: evidenceLevel
                evidence_type: evidenceType
                status
                variant_origin: variantOrigin
                molecular_profile: molecularProfile {
                  id
                }
                disease {
                  id
                }
                therapies {
                  id
                }
                phenotypes {
                  id
                }
                assertions {
                  id
                }
                source {
                    id
                }
                rating: evidenceRating
            }
        }"""



def _construct_get_all_evidence_payload():
    return """
        query evidenceItems($after: String) {
            evidence_items: evidenceItems(after: $after, status: ALL) {
                totalCount
                pageInfo {
                  hasNextPage
                  endCursor
                }
                nodes {
                    id
                    name
                    significance
                    description
                    therapy_interaction_type: therapyInteractionType
                    evidence_direction: evidenceDirection
                    evidence_level: evidenceLevel
                    evidence_type: evidenceType
                    status
                    variant_origin: variantOrigin
                    molecular_profile: molecularProfile {
                      id
                    }
                    disease {
                      id
                    }
                    therapies {
                      id
                    }
                    phenotypes {
                      id
                    }
                    assertions {
                      id
                    }
                    source {
                        id
                    }
                    rating: evidenceRating
                }
            }
        }"""



def _construct_get_assertion_payload():
    return """
        query assertion($id: Int!) {
            assertion(id: $id) {
                id
                name
                amp_level: ampLevel
                significance
                description
                therapy_interaction_type: therapyInteractionType
                assertion_direction: assertionDirection
                assertion_type: assertionType
                fda_companion_test: fdaCompanionTest
                fda_regulatory_approval: regulatoryApproval
                name
                nccn_guideline: nccnGuideline {
                  name
                }
                nccn_guideline_version: nccnGuidelineVersion
                status
                summary
                variant_origin: variantOrigin
                molecular_profile: molecularProfile {
                  id
                }
                acmg_codes: acmgCodes {
                  id
                  code
                  description
                }
                clingen_codes: clingenCodes {
                    id
                    code
                    description
                }
                disease {
                  id
                }
                endorsements {
                  nodes {
                    id
                  }
                }
                therapies {
                  id
                }
                evidenceItems {
                  id
                }
                phenotypes {
                  id
                }
            }
        }"""


def _construct_get_all_assertions_payload():
    return """
        query assertions($after: String) {
            assertions(after: $after, status: ALL) {
                totalCount
                pageInfo {
                  hasNextPage
                  endCursor
                }
                nodes {
                    id
                    name
                    amp_level: ampLevel
                    significance
                    description
                    therapy_interaction_type: therapyInteractionType
                    assertion_direction: assertionDirection
                    assertion_type: assertionType
                    fda_companion_test: fdaCompanionTest
                    fda_regulatory_approval: regulatoryApproval
                    name
                    nccn_guideline: nccnGuideline {
                      name
                    }
                    nccn_guideline_version: nccnGuidelineVersion
                    status
                    summary
                    variant_origin: variantOrigin
                    molecular_profile: molecularProfile {
                      id
                    }
                    acmg_codes: acmgCodes {
                      id
                      code
                      description
                    }
                    clingen_codes: clingenCodes {
                        id
                        code
                        description
                    }
                    disease {
                      id
                    }
                    endorsements {
                      nodes {
                        id
                      }
                    }
                    therapies {
                      id
                    }
                    evidenceItems {
                      id
                    }
                    phenotypes {
                      id
                    }
                }
            }
        }"""


def _construct_get_variant_group_payload():
    return """
        query variantGroup($id: Int!) {
            variant_group: variantGroup(id: $id) {
                id
                name
                description
                variants(first: 100) {
                  nodes {
                    id
                  }
                }
                sources {
                    id
                }
            }
        }"""


def _construct_get_all_variant_groups_payload():
    return """
        query variantGroups($after: String) {
            variant_groups: variantGroups(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                description
                variants(first: 100) {
                  nodes {
                    id
                  }
                }
                sources {
                    id
                }
              }
            }
        }"""


def _construct_get_source_payload():
    return """
        query source($id: Int!) {
            source(id: $id) {
                id
                name
                title
                citation
                citation_id: citationId
                source_type: sourceType
                abstract
                asco_abstract_id: ascoAbstractId
                author_string: authorString
                full_journal_title: fullJournalTitle
                journal
                pmc_id: pmcId
                publication_date: publicationDate
                source_url: sourceUrl
                clinical_trials: clinicalTrials {
                    id
                    name
                    description
                    nctId
                    url
                }
            }
        }"""


def _construct_get_all_sources_payload():
    return """
        query sources($after: String) {
            sources(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                title
                citation
                citation_id: citationId
                source_type: sourceType
                abstract
                asco_abstract_id: ascoAbstractId
                author_string: authorString
                full_journal_title: fullJournalTitle
                journal
                pmc_id: pmcId
                publication_date: publicationDate
                source_url: sourceUrl
                clinical_trials: clinicalTrials {
                    id
                    name
                    description
                    nctId
                    url
                }
              }
            }
        }"""


def _construct_get_disease_payload():
    return """
        query disease($id: Int!) {
            disease(id: $id) {
              id
              name
              doid
              disease_url: diseaseUrl
              aliases: diseaseAliases
            }
        }"""


def _construct_get_all_diseases_payload():
    return """
        query diseases($after: String) {
            diseases(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                  id
                  name
                  doid
                  disease_url: diseaseUrl
                  aliases: diseaseAliases
              }
            }
        }
    """


def _construct_get_therapy_payload():
    return """
        query therapy($id: Int!) {
            therapy(id: $id) {
              id
              name
              ncit_id: ncitId
              aliases: therapyAliases
              therapy_url: therapyUrl
            }
        }"""


def _construct_get_all_therapies_payload():
    return """
        query therapies($after: String) {
            therapies(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                  id
                  name
                  ncit_id: ncitId
                  aliases: therapyAliases
                  therapy_url: therapyUrl
              }
            }
        }
    """


def _construct_get_phenotype_payload():
    return """
        query phenotype($id: Int!) {
            phenotype(id: $id) {
              id
              name
              hpo_id: hpoId
              phenotype_url: url
            }
        }"""


def _construct_get_all_phenotypes_payload():
    return """
        query phenotypes($after: String) {
            phenotypes(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                  id
                  name
                  hpo_id: hpoId
                  phenotype_url: url
              }
            }
        }
    """

def _construct_get_organization_payload():
    return """
        query organization($id: Int!) {
            organization(id: $id) {
                id
                name
                url
                description
            }
        }"""


def _construct_get_all_organizations_payload():
    return """
        query organizations($after: String) {
            organizations(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                url
                description
              }
            }
        }
    """

def _construct_get_endorsement_payload():
    return """
        query endorsement($id: Int!) {
            endorsement(id: $id) {
                id
                assertion {
                  id
                }
                organization {
                  id
                }
                status
                last_reviewed: lastReviewed
                ready_for_clinvar_submission: readyForClinvarSubmission
            }
        }"""


def _construct_get_all_endorsements_payload():
    return """
        query endorsements($after: String) {
            endorsements(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                assertion {
                  id
                }
                organization {
                  id
                }
                status
                last_reviewed: lastReviewed
                ready_for_clinvar_submission: readyForClinvarSubmission
              }
            }
        }
    """
