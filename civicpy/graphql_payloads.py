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


def _construct_get_all_genes_payload():
    return """
        query genes($after: String) {
            genes(after: $after) {
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
                    feature {
                        id
                    }
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


def _construct_get_all_factors_payload():
    return """
        query factors($after: String) {
            factors(after: $after) {
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


def _construct_get_all_fusions_payload():
    return """
        query fusions($after: String) {
            fusions(after: $after) {
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
                    }
                    ... on Variant {
                        id
                        name
                        deprecated
                    }
                }
                sources {
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
                    }
                    ... on Variant {
                        id
                        name
                        deprecated
                    }
                }
                sources {
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
                  name
                  display_name: displayName
                  doid
                  disease_url: diseaseUrl
                  aliases: diseaseAliases
                }
                therapies {
                  id
                  name
                  ncit_id: ncitId
                  therapy_url: therapyUrl
                  aliases: therapyAliases
                }
                phenotypes {
                  id
                  name
                  hpo_id: hpoId
                  url
                }
                assertions {
                  id
                }
                source {
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
                      name
                      display_name: displayName
                      doid
                      disease_url: diseaseUrl
                      aliases: diseaseAliases
                    }
                    therapies {
                      id
                      name
                      ncit_id: ncitId
                      therapy_url: therapyUrl
                      aliases: therapyAliases
                    }
                    phenotypes {
                      id
                      name
                      hpo_id: hpoId
                      url
                    }
                    assertions {
                      id
                    }
                    source {
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
                  name
                  display_name: displayName
                  doid
                  disease_url: diseaseUrl
                  aliases: diseaseAliases
                }
                therapies {
                  id
                  name
                  ncit_id: ncitId
                  therapy_url: therapyUrl
                  aliases: therapyAliases
                }
                evidenceItems {
                  id
                }
                phenotypes {
                  id
                  name
                  hpo_id: hpoId
                  url
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
                      name
                      display_name: displayName
                      doid
                      disease_url: diseaseUrl
                      aliases: diseaseAliases
                    }
                    therapies {
                      id
                      name
                      ncit_id: ncitId
                      therapy_url: therapyUrl
                      aliases: therapyAliases
                    }
                    evidenceItems {
                      id
                    }
                    phenotypes {
                      id
                      name
                      hpo_id: hpoId
                      url
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
            }
        }"""
