import pytest
from civicpy import civic
from civicpy.civic import CoordinateQuery
import logging

ELEMENTS = [
    'assertion'
]


@pytest.fixture(scope="module", params=ELEMENTS)
def element(request):
    element_type = request.param
    return civic._get_elements_by_ids(element_type, [1])[0]


@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)

@pytest.fixture(scope="module")
def v600e_mp(v600e):
    return v600e.single_variant_molecular_profile


@pytest.fixture(scope="module")
def v600e_assertion():
    return civic.get_assertion_by_id(7)


class TestGetFunctions(object):
    @pytest.mark.skip(reason="Queries live API - Move to heartbeat test")
    def test_element_lookup_by_id(self):
        assertion = civic.element_lookup_by_id('assertion', '1')
        assert assertion['id'] == 1

    @pytest.mark.skip(reason="Queries live API - Move to heartbeat test")
    def test_get_assertions(self):
        test_ids = [1, 2, 3]
        results = civic._get_elements_by_ids('assertion', test_ids)
        assert len(results) == 3


class TestCivicRecord(object):

    def test_module(self):
        assert str(type(civic.MODULE)) == "<class 'module'>"


class TestElements(object):

    @pytest.mark.skip(reason="Queries live API - Move to heartbeat test")
    def test_attribute_fail(self, element):
        with pytest.raises(AttributeError):
            element.foo

    @pytest.mark.skip(reason="Queries live API - Move to heartbeat test")
    def test_completeness(self, element):
        for complex_field in element._COMPLEX_FIELDS:
            complex_value = getattr(element, complex_field)
            if not complex_value:
                continue
            if isinstance(complex_value, list):
                complex_value = complex_value[0]
            if isinstance(complex_value, civic.CivicAttribute):
                assert not complex_value._partial


class TestEvidence(object):

    def test_get_source_ids(self, v600e_mp):
        assert len(v600e_mp.evidence)
        assert len(v600e_mp.evidence) / 2 <= len(v600e_mp.evidence_sources)
        for source in v600e_mp.evidence_sources:
            assert source.citation_id
            assert source.source_type
            assert hasattr(source, 'abstract')
            assert hasattr(source, 'asco_abstract_id')
            assert source.author_string != None
            assert source.citation
            assert source.full_journal_title
            assert source.journal
            assert hasattr(source, 'pmc_id')
            assert source.publication_date
            assert source.source_url
            assert source.title
            assert hasattr(source, 'clinical_trials')

    def test_get_all(self):
        evidence = civic.get_all_evidence()
        assert len(evidence) >= 6656

    def test_get_by_id(self):
        evidence = civic.get_evidence_by_id(1637)
        assert evidence.type == 'evidence'
        assert evidence.id == 1637

    def test_get_non_rejected(self):
        evidence = civic.get_all_evidence(include_status=['accepted', 'submitted'])
        assert len(evidence) >= 6519

    def test_get_accepted_only(self):
        evidence = civic.get_all_evidence(include_status=['accepted'])
        assert len(evidence) >= 3247

    def test_properties(self, v600e_mp):
        evidence = v600e_mp.evidence[0]
        assert evidence.molecular_profile.name == 'BRAF V600E'
        assert evidence.statement == evidence.description
        assert len(evidence.assertions) == 0
        assert evidence.disease.name == 'Skin Melanoma'
        assert len(evidence.therapies) == 1
        assert len(evidence.phenotypes) == 0

        evidence = civic.get_evidence_by_id(7365)
        expected_phenotype_ids = {15320, 16642}
        assert set(evidence.phenotype_ids) == expected_phenotype_ids
        assert len(evidence.phenotypes) == 2
        for p in evidence.phenotypes:
            assert p.id in expected_phenotype_ids

class TestVariants(object):

    def test_get_all(self):
        variants = civic.get_all_variants()
        assert len(variants) >= 3815

    def test_get_non_rejected(self):
        variants = civic.get_all_variants(include_status=['accepted', 'submitted'])
        assert len(variants) >= 2374

    def test_get_accepted_only(self):
        variants = civic.get_all_variants(include_status=['accepted'])
        assert len(variants) >= 1333

    def test_get_by_id(self):
        variant = civic.get_variant_by_id(12)
        assert variant.type == 'variant'
        assert variant.id == 12

    def test_get_by_name(self, v600e):
        variants = civic.search_variants_by_name("V600E")
        assert len(variants) == 1
        assert variants[0] == v600e

    def test_get_by_caid(self, v600e):
        variants = civic.search_variants_by_allele_registry_id("CA123643")
        assert len(variants) == 1
        assert variants[0] == v600e

    def test_get_by_hgvs(self, v600e):
        variants = civic.search_variants_by_hgvs("ENST00000288602.6:c.1799T>A")
        assert len(variants) == 1
        assert variants[0] == v600e

    def test_sanitize_coordinate_bases(self):
        variant1 = civic.get_variant_by_id(2696)
        variant2 = civic.get_variant_by_id(558)
        for v in variant1, variant2:
            assert v.coordinates.reference_bases not in ['', '-']
            assert v.coordinates.variant_bases not in ['', '-']

    def test_shared_properties(self):
        variant = civic.get_variant_by_id(11)
        assert sorted(variant.aliases) == sorted(variant.variant_aliases)
        assert sorted(variant.groups) == sorted(variant.variant_groups)
        assert sorted(variant.types) == sorted(variant.variant_types)
        assert len(variant.molecular_profiles) == 1
        assert variant.single_variant_molecular_profile.id == 11


class TestGeneVariants(object):

    def test_get_all(self):
        variants = civic.get_all_gene_variants()
        assert len(variants) >= 3522
        for variant in variants:
            assert variant.subtype == 'gene_variant'

    def test_get_by_id(self):
        variant = civic.get_variant_by_id(11)
        assert variant.id == 11
        assert variant.type == 'variant'

    def test_attributes(self):
        variant = civic.get_variant_by_id(11)
        assert variant.coordinates.ensembl_version == 75
        assert variant.entrez_name == "BRAF"
        assert variant.entrez_id == 673

    def test_properties(self):
        variant = civic.get_variant_by_id(11)
        assert variant.gene.id == 5
        assert variant.gene == variant.feature
        assert variant.is_insertion == False
        assert variant.is_deletion == False


class TestFusionVariants(object):

    def test_get_all(self):
        variants = civic.get_all_fusion_variants()
        assert len(variants) >= 263
        for variant in variants:
            assert variant.subtype == 'fusion_variant'

    def test_get_by_id(self):
        variant = civic.get_variant_by_id(1)
        assert variant.id == 1
        assert variant.type == 'variant'

    def test_attributes(self):
        variant = civic.get_variant_by_id(1)
        assert variant.vicc_compliant_name == 'BCR(entrez:613)::ABL1(entrez:25)'
        assert variant.five_prime_start_exon_coordinates.chromosome == "22"
        assert variant.five_prime_start_exon_coordinates.reference_build == "GRCH37"
        assert variant.five_prime_start_exon_coordinates.ensembl_id == "ENSE00001897802"
        assert variant.five_prime_start_exon_coordinates.ensembl_version == 75
        assert variant.five_prime_start_exon_coordinates.representative_transcript == "ENST00000305877.8"
        assert variant.five_prime_start_exon_coordinates.strand == "POSITIVE"
        assert variant.five_prime_start_exon_coordinates.exon == 1
        assert variant.five_prime_start_exon_coordinates.exon_offset==0
        assert variant.five_prime_start_exon_coordinates.exon_offset_direction is None
        assert variant.five_prime_start_exon_coordinates.start == 23522397
        assert variant.five_prime_start_exon_coordinates.stop == 23524426

        assert variant.five_prime_end_exon_coordinates.chromosome == "22"
        assert variant.five_prime_end_exon_coordinates.reference_build == "GRCH37"
        assert variant.five_prime_end_exon_coordinates.ensembl_id == "ENSE00001781765"
        assert variant.five_prime_end_exon_coordinates.ensembl_version == 75
        assert variant.five_prime_end_exon_coordinates.representative_transcript == "ENST00000305877.8"
        assert variant.five_prime_end_exon_coordinates.strand == "POSITIVE"
        assert variant.five_prime_end_exon_coordinates.exon == 14
        assert variant.five_prime_end_exon_coordinates.exon_offset==0
        assert variant.five_prime_end_exon_coordinates.exon_offset_direction is None
        assert variant.five_prime_end_exon_coordinates.start == 23632526
        assert variant.five_prime_end_exon_coordinates.stop == 23632600

        assert variant.three_prime_start_exon_coordinates.chromosome == "9"
        assert variant.three_prime_start_exon_coordinates.reference_build == "GRCH37"
        assert variant.three_prime_start_exon_coordinates.ensembl_id == "ENSE00000984287"
        assert variant.three_prime_start_exon_coordinates.ensembl_version == 75
        assert variant.three_prime_start_exon_coordinates.representative_transcript == "ENST00000318560.5"
        assert variant.three_prime_start_exon_coordinates.strand == "POSITIVE"
        assert variant.three_prime_start_exon_coordinates.exon == 2
        assert variant.three_prime_start_exon_coordinates.exon_offset==0
        assert variant.three_prime_start_exon_coordinates.exon_offset_direction is None
        assert variant.three_prime_start_exon_coordinates.start == 133729451
        assert variant.three_prime_start_exon_coordinates.stop == 133729624

        assert variant.three_prime_end_exon_coordinates.chromosome == "9"
        assert variant.three_prime_end_exon_coordinates.reference_build == "GRCH37"
        assert variant.three_prime_end_exon_coordinates.ensembl_id == "ENSE00001457584"
        assert variant.three_prime_end_exon_coordinates.ensembl_version == 75
        assert variant.three_prime_end_exon_coordinates.representative_transcript == "ENST00000318560.5"
        assert variant.three_prime_end_exon_coordinates.strand == "POSITIVE"
        assert variant.three_prime_end_exon_coordinates.exon == 11
        assert variant.three_prime_end_exon_coordinates.exon_offset==0
        assert variant.three_prime_end_exon_coordinates.exon_offset_direction is None
        assert variant.three_prime_end_exon_coordinates.start == 133759356
        assert variant.three_prime_end_exon_coordinates.stop == 133763062

    def test_nullable_fields(self):
        variant = civic.get_variant_by_id(5041)
        assert variant.three_prime_coordinates is None
        assert variant.three_prime_start_exon_coordinates is None
        assert variant.three_prime_end_exon_coordinates is None

    def test_properties(self):
        variant = civic.get_variant_by_id(1)
        assert variant.fusion.id == 62007
        assert variant.fusion == variant.feature


class TestFactorVariants(object):

    def test_get_all(self):
        variants = civic.get_all_factor_variants()
        assert len(variants) >= 8
        for variant in variants:
            assert variant.subtype == 'factor_variant'

    def test_get_by_id(self):
        variant = civic.get_variant_by_id(4985)
        assert variant.id == 4985
        assert variant.type == 'variant'

    def test_attributes(self):
        variant = civic.get_variant_by_id(4985)
        assert variant.ncit_id == 'C131459'

    def test_properties(self):
        variant = civic.get_variant_by_id(4985)
        assert variant.factor.id == 61746
        assert variant.factor == variant.feature


class TestMolecularProfiles(object):

    def test_get_all(self):
        mps = civic.get_all_molecular_profiles()
        assert len(mps) >= 2396

    def test_get_non_rejected(self):
        mps = civic.get_all_molecular_profiles(include_status=['accepted', 'submitted'])
        assert len(mps) >= 2374

    def test_get_accepted_only(self):
        mps = civic.get_all_molecular_profiles(include_status=['accepted'])
        assert len(mps) >= 1333

    def test_get_by_id(self):
        mp = civic.get_molecular_profile_by_id(12)
        assert mp.type == 'molecular_profile'
        assert mp.id == 12

    def test_get_by_id_complex_mp(self):
        mp = civic.get_molecular_profile_by_id(4432)
        assert mp.type == 'molecular_profile'
        mp_parsed_name = mp.parsed_name
        assert len(mp_parsed_name) == 5
        egfr_gene = mp_parsed_name[0]
        assert egfr_gene.type == "gene"
        assert egfr_gene.id == 19
        assert egfr_gene.name == "EGFR"
        variant0 = mp_parsed_name[1]
        assert variant0.type == "variant"
        assert variant0.id == 33
        assert variant0.name == "L858R"
        text_segment = mp_parsed_name[2]
        assert text_segment == "OR"
        assert mp_parsed_name[3] == egfr_gene
        variant1 = mp_parsed_name[4]
        assert variant1.type == "variant"
        assert variant1.id == 133
        assert variant1.name == "Exon 19 Deletion"

    def test_properties(self):
        mp = civic.get_molecular_profile_by_id(4432)
        assert len(mp.evidence_sources) == 10
        assert mp.summary == mp.description
        assert len(mp.evidence_items) == 11
        assert len(mp.assertions) == 0
        assert len(mp.variants) == 2
        assert len(mp.sources) == 0


class TestVariantGroups(object):

    def test_get_all(self):
        variant_groups = civic.get_all_variant_groups()
        assert len(variant_groups) >= 25

    def test_get_by_id(self):
        variant_group = civic.get_variant_group_by_id(1)
        assert variant_group.type == 'variant_group'
        assert variant_group.id == 1

    def test_properties(self):
        variant_group = civic.get_variant_group_by_id(1)
        assert len(variant_group.variants) == 7
        assert len(variant_group.sources) == 0


class TestAssertions(object):

    def test_get_all(self):
        assertions = civic.get_all_assertions()
        assert len(assertions) >= 28

    def test_get_non_rejected(self):
        assertions = civic.get_all_assertions(include_status=['accepted', 'submitted'])
        assert len(assertions) >= 24

    def test_get_accepted_only(self):
        assertions = civic.get_all_assertions(include_status=['accepted'])
        assert len(assertions) >= 16

    def test_get_by_id(self):
        assertion = civic.get_assertion_by_id(7)
        assert assertion.type == 'assertion'
        assert assertion.id == 7

    def test_has_boolean_flags(self, v600e_assertion):
        assert v600e_assertion.fda_companion_test is True
        assert v600e_assertion.fda_regulatory_approval is True

    def test_properties(self):
        assertion = civic.get_assertion_by_id(18)
        assert assertion.evidence == assertion.evidence_items
        assert assertion.hpo_ids == [p.hpo_id for p in assertion.phenotypes if p.hpo_id]
        assert assertion.acmg_codes
        for acmg_code in assertion.acmg_codes:
            assert acmg_code.id
            assert acmg_code.code
            assert acmg_code.description
        assert assertion.disease.name == "Von Hippel-Lindau Disease"
        assert len(assertion.therapies) == 0
        assert len(assertion.phenotypes) == 3
        assert len(assertion.endorsements) == 0
        assert assertion.molecular_profile.id == 1686

        # Test assertion with clingen_codes
        assertion = civic.get_assertion_by_id(53)
        assert assertion.clingen_codes
        for clingen_code in assertion.clingen_codes:
            assert clingen_code.id
            assert clingen_code.code
            assert clingen_code.description


class TestFeatures(object):

    def test_get_all(self):
        features = civic.get_all_features()
        assert len(features) >= 407

    def test_get_non_rejected(self):
        features = civic.get_all_features(include_status=['accepted', 'submitted'])
        assert len(features) >= 402

    def test_get_accepted_only(self):
        features = civic.get_all_features(include_status=['accepted'])
        assert len(features) >= 322

    def test_get_by_id(self):
        feature = civic.get_feature_by_id(58)
        assert feature.type == 'gene'
        assert feature.id == 58

    def test_get_by_ids(self):
        features = civic.get_features_by_ids([58, 61748, 61753])
        assert features[0].type == 'gene'
        assert features[0].id == 58
        assert features[1].type == 'factor'
        assert features[1].id == 61748
        assert features[2].type == 'fusion'
        assert features[2].id == 61753

    def test_attributes(self):
        feature = civic.get_feature_by_id(58)
        assert feature.name == 'VHL'


class TestGenes(object):

    def test_get_all(self):
        genes = civic.get_all_genes()
        assert len(genes) >= 407

    def test_get_non_rejected(self):
        genes = civic.get_all_genes(include_status=['accepted', 'submitted'])
        assert len(genes) >= 402

    def test_get_accepted_only(self):
        genes = civic.get_all_genes(include_status=['accepted'])
        assert len(genes) >= 322

    def test_get_by_id(self):
        gene = civic.get_gene_by_id(58)
        assert gene.type == 'gene'
        assert gene.id == 58

    def test_get_by_name(self):
        gene = civic.get_gene_by_name('BRAF')
        assert gene.name == 'BRAF'
        assert gene.type == 'gene'

    def test_get_by_entrez_id(self):
        gene = civic.get_gene_by_entrez_id(673)
        assert gene.entrez_id == 673
        assert gene.type == 'gene'

    def test_attributes(self):
        gene = civic.get_gene_by_id(58)
        assert gene.name == 'VHL'

    def test_properties(self):
        gene = civic.get_gene_by_id(58)
        assert len(gene.variants) == 862
        assert len(gene.sources) == 4


class TestFactors(object):

    def test_get_all(self):
        factors = civic.get_all_factors()
        assert len(factors) >= 6

    def test_get_non_rejected(self):
        factors = civic.get_all_factors(include_status=['accepted', 'submitted'])
        assert len(factors) >= 6

    def test_get_accepted_only(self):
        factors = civic.get_all_factors(include_status=['accepted'])
        assert len(factors) >= 2

    def test_get_by_id(self):
        factor = civic.get_factor_by_id(61748)
        assert factor.type == 'factor'
        assert factor.id == 61748

    def test_get_by_name(self):
        factor = civic.get_factor_by_name('MSI')
        assert factor.type == 'factor'
        assert factor.name == 'MSI'

    def test_get_by_ncit_id(self):
        factor = civic.get_factor_by_ncit_id('C36318')
        assert factor.type == 'factor'
        assert factor.ncit_id == 'C36318'

    def test_attributes(self):
        factor = civic.get_factor_by_id(61748)
        assert factor.name == 'CK'
        assert factor.full_name == 'Complex Karyotype'

    def test_properties(self):
        factor = civic.get_factor_by_id(61748)
        assert len(factor.variants) == 1
        assert len(factor.sources) == 0


class TestFusions(object):

    def test_get_all(self):
        fusions = civic.get_all_fusions()
        assert len(fusions) >= 256

    def test_get_non_rejected(self):
        fusions = civic.get_all_fusions(include_status=['accepted', 'submitted'])
        assert len(fusions) >= 255

    def test_get_accepted_only(self):
        fusions = civic.get_all_fusions(include_status=['accepted'])
        assert len(fusions) >= 166

    def test_get_by_id(self):
        fusion = civic.get_fusion_by_id(61753)
        assert fusion.type == 'fusion'
        assert fusion.id == 61753

    def test_attributes(self):
        fusion = civic.get_fusion_by_id(61753)
        assert fusion.name == 'MEF2D::CSF1R'
        assert fusion.five_prime_gene.name == 'MEF2D'
        assert fusion.three_prime_gene.name == 'CSF1R'

    def test_properties(self):
        fusion = civic.get_fusion_by_id(61753)
        assert len(fusion.variants) == 1
        assert len(fusion.sources) == 0

    def test_get_fusion_by_name(self):
        fusion = civic.get_fusion_by_name("BCR::ABL1")
        assert fusion.id == 62007
        assert fusion.name == 'BCR::ABL1'

    def test_search_fusions_by_partner_gene_id(self):
        fusions = civic.search_fusions_by_partner_gene_id(573)
        assert len(fusions) >= 5


class TestDiseases(object):

    def test_get_all(self):
        diseases = civic.get_all_diseases()
        assert len(diseases) >= 418

    def test_get_by_id(self):
        d = civic.get_disease_by_id(22)
        assert d.id == 22
        assert d.type == 'disease'

    def test_attributes(self):
        breast_cancer = civic.get_disease_by_id(22)
        assert breast_cancer.doid == '1612'
        assert breast_cancer.name == 'Breast Cancer'
        assert set(breast_cancer.aliases) == {
            'Breast Tumor',
            'Malignant Neoplasm Of Breast',
            'Malignant Tumor Of The Breast',
            'Mammary Cancer',
            'Mammary Tumor',
            'Primary Breast Cancer'
        }

    def test_get_by_name(self):
        breast_cancer = civic.get_disease_by_name('Breast Cancer')
        assert breast_cancer.id == 22

    def test_get_by_doid(self):
        breast_cancer = civic.get_disease_by_doid('1612')
        assert breast_cancer.id == 22

    def test_properties(self):
        breast_cancer = civic.get_disease_by_id(22)
        assert len(breast_cancer.evidence) >= 278
        assert breast_cancer.evidence == breast_cancer.evidence_items
        assert len(breast_cancer.assertions) >= 2


class TestTherapies(object):

    def test_get_all(self):
        therapies = civic.get_all_therapies()
        assert len(therapies) >= 552

    def test_get_by_id(self):
        t = civic.get_therapy_by_id(19)
        assert t.id == 19
        assert t.type == 'therapy'

    def test_attributes(self):
        trametinib = civic.get_therapy_by_id(19)
        assert trametinib.ncit_id == 'C77908'
        assert trametinib.name == 'Trametinib'
        assert set(trametinib.aliases) == {
            'JTP-74057',
            'GSK1120212',
            'MEK Inhibitor GSK1120212',
            'Mekinist',
            'N-(3-{3-cyclopropyl-5-[(2-fluoro-4-iodophenyl)amino]-6,8-dimethyl-2,4,7-trioxo-3,4,6,7-tetrahydropyrido[4,3-d]pyrimidin-1(2H)-yl}phenyl)acetamide'
        }

    def test_get_by_name(self):
        trametinib = civic.get_therapy_by_name('Trametinib')
        assert trametinib.id == 19

    def test_get_by_ncit_id(self):
        trametinib = civic.get_therapy_by_ncit_id('C77908')
        assert trametinib.id == 19

    def test_properties(self):
        trametinib = civic.get_therapy_by_id(19)
        assert len(trametinib.evidence) >= 136
        assert trametinib.evidence == trametinib.evidence_items
        assert len(trametinib.assertions) >= 3


class TestSource(object):
    def test_get_all(self):
        sources = civic.get_all_sources()
        assert len(sources) >= 3868

    def test_get_by_id(self):
        s = civic.get_source_by_id(1)
        assert s.id == 1
        assert s.type == 'source'

    def test_attributes(self):
        s = civic.get_source_by_id(947)
        assert s.citation == 'McArthur et al., 2014'
        assert s.citation_id == '24508103'
        assert s.source_type == 'PUBMED'
        assert s.abstract.startswith('In the BRIM-3 trial,')
        assert s.author_string.startswith('Grant A McArthur,')
        assert s.full_journal_title == 'The Lancet. Oncology'
        assert s.journal == 'Lancet Oncol'
        assert s.pmc_id == 'PMC4382632'
        assert s.publication_date == '2014-3'
        assert s.source_url == 'http://www.ncbi.nlm.nih.gov/pubmed/24508103'
        assert s.title.startswith('Safety and efficacy of vemurafenib')
        assert len(s.clinical_trials) == 1

    def test_get_pubmed_source_by_id(self):
        s = civic.get_pubmed_source_by_id('24889366')
        assert s.citation_id == '24889366'
        assert s.source_type == 'PUBMED'

    def test_get_ash_source_by_doi(self):
        s = civic.get_ash_source_by_doi('10.1182/blood-2021-145491')
        assert s.citation_id == '10.1182/blood-2021-145491'
        assert s.source_type == 'ASH'

    def test_get_asco_source_by_id(self):
        s = civic.get_asco_source_by_id('144555')
        assert s.citation_id == '144555'
        assert s.asco_abstract_id == 5005
        assert s.source_type == 'ASCO'


class TestOrganization(object):
    def test_get_all(self):
        organizations = civic.get_all_organizations()
        assert len(organizations) >= 23

    def test_get_by_id(self):
        o = civic.get_organization_by_id(1)
        assert o.id == 1
        assert o.type == 'organization'

    def test_attributes(self):
        org = civic.get_organization_by_id(1)
        assert org.name == 'The McDonnell Genome Institute'
        assert org.url == 'http://genome.wustl.edu/'

    def test_properties(self):
        org = civic.get_organization_by_id(14)
        assert len(org.endorsements) >= 1


class TestEndorsement(object):
    def test_get_all(self):
        endorsements = civic.get_all_endorsements()
        assert len(endorsements) >= 4

    def test_get_by_id(self):
        e = civic.get_endorsement_by_id(1)
        assert e.id == 1
        assert e.type == 'endorsement'

    def test_attributes(self):
        e = civic.get_endorsement_by_id(1)
        assert e.organization_id == 14
        assert e.assertion_id == 101

    def test_search_endorsements_by_organization_id(self):
        endorsements = civic.search_endorsements_by_organization_id(1)
        assert len(endorsements) >= 1

    def test_search_endorsements_by_assertion_id(self):
        endorsements = civic.search_endorsements_by_assertion_id(6)
        assert len(endorsements) >= 1


class TestPhenotypes(object):

    def test_get_all(self):
        phenotypes = civic.get_all_phenotypes()
        assert len(phenotypes) >= 265

    def test_get_by_id(self):
        pediatric_onset = civic.get_phenotype_by_id(15320)
        assert pediatric_onset.hpo_id == 'HP:0410280'
        assert pediatric_onset.name == 'Pediatric onset'

    def test_get_by_name(self):
        pediatric_onset = civic.get_phenotype_by_name('Pediatric onset')
        assert pediatric_onset.id == 15320

    def test_get_by_hpo_id(self):
        pediatric_onset = civic.get_phenotype_by_hpo_id('HP:0410280')
        assert pediatric_onset.id == 15320

    def test_properties(self):
        pediatric_onset = civic.get_phenotype_by_id(15320)
        assert len(pediatric_onset.evidence) >= 128
        assert pediatric_onset.evidence == pediatric_onset.evidence_items
        assert len(pediatric_onset.assertions) >= 25


class TestCoordinateSearch(object):

    def test_search_assertions(self):
        query = CoordinateQuery('7', 140453136, 140453136, 'T', '*')
        assertions = civic.search_assertions_by_coordinates(query)
        assertion_ids = [x.id for x in assertions]
        v600e_assertion_ids = (7, 10, 12, 20, 23)
        v600k_assertion_ids = (11, 13)
        assert set(assertion_ids) >= set(v600e_assertion_ids + v600k_assertion_ids)
        assertions = civic.search_assertions_by_coordinates(query, search_mode='exact')
        assertion_ids = [x.id for x in assertions]
        assert set(assertion_ids) >= set(v600e_assertion_ids)

    def test_search_evidence(self):
        pass

    def test_single_and_bulk_exact_return_same_variants(self):
        query = CoordinateQuery('7', 140453136, 140453136, 'T', '*')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('7', 140453136, 140453136, 'T', 'A')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('7', 140453136, 140453136, 'T', None)
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 0
        assert len(variants_bulk) == 0

        query = CoordinateQuery('7', 140453136, 140453137, 'TT', '*')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('7', 140453136, 140453137, 'TT', 'AC')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('7', 140453136, 140453137, 'TT', None)
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 0
        assert len(variants_bulk) == 0

        query = CoordinateQuery('3', 10183694, 10183694, None, 'G')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('3', 10183694, 10183694, 'T', 'G')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('3', 10183694, 10183694, '*', 'G')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        variants_single = list(map(lambda v: hash(v), variants_single))
        variants_bulk = list(map(lambda v: v.v_hash, variants_bulk[query]))
        assert len(variants_single) == 2
        assert len(variants_bulk) == 2
        assert sorted(variants_single) == sorted(variants_bulk)
        assert sorted(variants_single) == sorted(variants_bulk)

    def test_bulk_any_search_variants(self):
        sorted_queries = [
            CoordinateQuery('7', 140453136, 140453136, 'T'),
            CoordinateQuery('7', 140453136, 140453137, 'TT')
        ]
        search_results = civic.bulk_search_variants_by_coordinates(sorted_queries, search_mode='any')
        assert len(search_results[sorted_queries[0]]) >= 17
        assert len(search_results[sorted_queries[1]]) >= 11

    def test_bulk_exact_search_variants(self):
        sorted_queries = [
            CoordinateQuery('7', 140453136, 140453136, 'T', '*'),
            CoordinateQuery('7', 140453136, 140453137, 'TT', '*'),
            CoordinateQuery('7', 140453136, 140453136, 'T', 'A'),
            CoordinateQuery('7', 140453136, 140453137, 'TT', 'AC'),
        ]
        search_results = civic.bulk_search_variants_by_coordinates(sorted_queries, search_mode='exact')
        assert len(search_results[sorted_queries[0]]) == 1
        assert len(search_results[sorted_queries[1]]) == 1
        assert len(search_results[sorted_queries[2]]) == 1
        assert len(search_results[sorted_queries[3]]) == 1

    def test_bulk_qe_search_variants(self):
        sorted_queries = [
            CoordinateQuery('7', 140453136, 140453136),
            CoordinateQuery('7', 140453136, 140453137)
        ]
        search_results = civic.bulk_search_variants_by_coordinates(sorted_queries, search_mode='query_encompassing')
        assert len(search_results[sorted_queries[0]]) == 1
        assert len(search_results[sorted_queries[1]]) >= 5

    def test_bulk_re_search_variants(self):
        sorted_queries = [
            CoordinateQuery('7', 140453136, 140453136),
            CoordinateQuery('7', 140453136, 140453137)
        ]
        search_results = civic.bulk_search_variants_by_coordinates(sorted_queries, search_mode='record_encompassing')
        assert len(search_results[sorted_queries[0]]) >= 17
        assert len(search_results[sorted_queries[1]]) >= 14

    def test_build38_exact_search_variants(self, v600e):
        query = CoordinateQuery('7', 140753336, 140753336, 'T', 'A', 'GRCh38')
        search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert len(search_results) == 1
        assert search_results[0] == v600e

        query = CoordinateQuery('7', 140753336, 140753337, 'TT', 'AC', 'GRCh38')
        search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert len(search_results) == 1
        assert search_results[0].id == 563

        query = CoordinateQuery('3', 10146548, 10146549, 'C', None, 'GRCh38')
        search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert len(search_results) >= 1
        assert search_results[0].id == 1918

        query = CoordinateQuery('3', 10146618, 10146618, None, 'G', 'GRCh38')
        search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert len(search_results) == 1
        assert search_results[0].id == 2042

    def test_build36_exact_search_variants(self, v600e):
        query = CoordinateQuery('7', 140099605, 140099605, 'T', 'A', 'NCBI36')
        search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert len(search_results) == 1
        assert search_results[0] == v600e

    def test_errors(self):
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, 'T', 'A')
            variants_single = civic.search_variants_by_coordinates(query, search_mode='wrong_mode')
        assert "unexpected search mode" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140753336, 140753336, '*', 'A', 'GRCh38')
            search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert "Can't use wildcard when searching for non-GRCh37 coordinates" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140753336, 140753336, None, None, 'GRCh38')
            search_results = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert "alt or ref required for non-GRCh37 coordinate queries" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140753336, 140753336, 'T', 'A', 'GRCh38')
            search_results = civic.search_variants_by_coordinates(query, search_mode='any')
        assert "Only exact search mode is supported for non-GRCh37 coordinate queries" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, '-', 'A')
            variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert "Unexpected alt `-` in coordinate query. Did you mean `None`?" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, 'T', '-')
            variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert "Unexpected ref `-` in coordinate query. Did you mean `None`?" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, '-', 'A', 'GRCh38')
            variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert "Unexpected alt `-` in coordinate query. Did you mean `None`?" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, 'T', '-', 'GRCh38')
            variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        assert "Unexpected ref `-` in coordinate query. Did you mean `None`?" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, '-', 'A')
            variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert "Unexpected alt `-` in coordinate query. Did you mean `None`?" in str(context.value)
        with pytest.raises(ValueError) as context:
            query = CoordinateQuery('7', 140453136, 140453136, 'T', '-')
            variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert "Unexpected ref `-` in coordinate query. Did you mean `None`?" in str(context.value)


#warning logging tests
LOGGER = logging.getLogger(__name__)

def test_is_valid_for_vcf_warnings(caplog):
    incomplete_coordinates_variant = civic.get_variant_by_id(27)
    incomplete_coordinates_variant.is_valid_for_vcf(emit_warnings=True)
    assert "Incomplete coordinates for variant 27. Skipping." in caplog.text

    unsupported_var_bases_variant = civic.get_variant_by_id(613)
    unsupported_var_bases_variant.is_valid_for_vcf(emit_warnings=True)
    assert "Unsupported variant base(s) for variant 613. Skipping." in caplog.text

    #currently no case for unsupported ref bases

def test_is_valid_for_gks_warnings_assertion(caplog):
    """Test that is_valid_for_gks_json works correctly for assertions"""
    not_accepted = civic.get_assertion_by_id(117)
    assert not not_accepted.is_valid_for_gks_json(emit_warnings=True)
    assert "Assertion 117 does not have 'accepted' status. Skipping" in caplog.text

    oncogenic_fusion = civic.get_assertion_by_id(101)
    assert not oncogenic_fusion.is_valid_for_gks_json(emit_warnings=True)
    assert "Assertion 101 type is not one of: 'DIAGNOSTIC', 'PREDICTIVE', or 'PROGNOSTIC'. Skipping" in caplog.text
    assert "Assertion 101 variant is not a ``GeneVariant``. Skipping" in caplog.text

    complex_mp = civic.get_assertion_by_id(88)
    assert not complex_mp.is_valid_for_gks_json(emit_warnings=True)
    assert "Assertion 88 has a complex molecular profile. Skipping" in caplog.text

def test_is_valid_for_gks_warnings_evidence(caplog):
    """Test that is_valid_for_gks_json works correctly for evidence items"""
    not_accepted_oncogenic_fusion = civic.get_evidence_by_id(6936)
    assert not not_accepted_oncogenic_fusion.is_valid_for_gks_json(emit_warnings=True)
    assert "Evidence 6936 type is not one of: 'DIAGNOSTIC', 'PREDICTIVE', or 'PROGNOSTIC'. Skipping" in caplog.text
    assert "Evidence 6936 variant is not a ``GeneVariant``. Skipping" in caplog.text
    assert "Evidence 6936 does not have 'accepted' status. Skipping" in caplog.text

    complex_mp = civic.get_evidence_by_id(8177)
    assert not complex_mp.is_valid_for_gks_json(emit_warnings=True)
    assert "Evidence 8177 has a complex molecular profile. Skipping" in caplog.text
