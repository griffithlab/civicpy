import pytest
from civicpy import civic, TEST_CACHE_PATH
from civicpy.civic import CoordinateQuery
import logging

ELEMENTS = [
    'Assertion'
]


def setup_module():
    civic.load_cache(local_cache_path=TEST_CACHE_PATH, on_stale='ignore')


@pytest.fixture(scope="module", params=ELEMENTS)
def element(request):
    element_type = request.param
    return civic._get_elements_by_ids(element_type, [1])[0]


@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)


@pytest.fixture(scope="module")
def v600e_assertion():
    return civic.get_assertion_by_id(7)


class TestGetFunctions(object):
    def test_element_lookup_by_id(self):
        assertion = civic.element_lookup_by_id('assertion', '1')
        assert assertion['id'] == 1

    def test_get_assertions(self):
        test_ids = [1, 2, 3]
        results = civic._get_elements_by_ids('assertion', test_ids)
        assert len(results) == 3


class TestCivicRecord(object):

    def test_module(self):
        assert str(type(civic.MODULE)) == "<class 'module'>"


class TestElements(object):

    def test_attribute_fail(self, element):
        with pytest.raises(AttributeError):
            element.foo

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

    def test_get_source_ids(self, v600e):
        assert len(v600e.evidence)
        assert len(v600e.evidence) / 2 <= len(v600e.evidence_sources)
        for source in v600e.evidence_sources:
            assert source.citation_id
            assert source.source_type

    def test_get_all(self):
        evidence = civic.get_all_evidence()
        assert len(evidence) >= 6656

    def test_get_non_rejected(self):
        evidence = civic.get_all_evidence(include_status=['accepted', 'submitted'])
        assert len(evidence) >= 6519

    def test_get_accepted_only(self):
        evidence = civic.get_all_evidence(include_status=['accepted'])
        assert len(evidence) >= 3247

    #get_all_ids pulls from the live site so it will return more results than get_all_x
    def test_get_all_ids(self):
        evidence_ids = civic.get_all_evidence_ids()
        assert len(evidence_ids) >= len(civic.get_all_evidence())

    def test_properties(self, v600e):
        evidence = v600e.evidence[0]
        assert evidence.variant.name == 'V600E'
        assert evidence.statement == evidence.description

class TestVariants(object):

    def test_get_all(self):
        variants = civic.get_all_variants()
        assert len(variants) >= 2396

    def test_get_non_rejected(self):
        variants = civic.get_all_variants(include_status=['accepted', 'submitted'])
        assert len(variants) >= 2374

    def test_get_accepted_only(self):
        variants = civic.get_all_variants(include_status=['accepted'])
        assert len(variants) >= 1333

    #get_all_ids pulls from the live site so it will return more results than get_all_x
    def test_get_all_ids(self):
        variant_ids = civic.get_all_variant_ids()
        assert len(variant_ids) >= len(civic.get_all_variants())

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

    def test_properties(self):
        variant = civic.get_variant_by_id(11)
        assert sorted(variant.aliases) == sorted(variant.variant_aliases)
        assert sorted(variant.groups) == sorted(variant.variant_groups)
        assert sorted(variant.types) == sorted(variant.variant_types)
        assert variant.summary == variant.description


class TestVariantGroups(object):

    def test_get_all(self):
        variant_groups = civic.get_all_variant_groups()
        assert len(variant_groups) >= 25


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

    def test_has_boolean_flags(self, v600e_assertion):
        assert v600e_assertion.fda_companion_test is True
        assert v600e_assertion.fda_regulatory_approval is True

    def test_properties(self):
        assertion = civic.get_assertion_by_id(18)
        assert assertion.evidence == assertion.evidence_items
        assert assertion.hpo_ids == [p.hpo_id for p in assertion.phenotypes if p.hpo_id]


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

    #get_all_ids pulls from the live site so it might return more results than get_all_x
    def test_get_all_ids(self):
        genes = civic.get_all_genes(include_status=['accepted'])
        gene_ids = civic.get_all_gene_ids()
        assert len(gene_ids) >= len(genes)

    def test_get_by_id(self):
        gene = civic.get_gene_by_id(58)
        assert gene.name == 'VHL'

class TestCoordinateSearch(object):

    def test_search_assertions(self):
        query = CoordinateQuery('7', 140453136, 140453136, 'T', '*')
        assertions = civic.search_assertions_by_coordinates(query)
        assertion_ids = [x.id for x in assertions]
        v600e_assertion_ids = (7, 10, 12, 20)
        v600k_assertion_ids = (11, 13)
        assert set(assertion_ids) >= set(v600e_assertion_ids + v600k_assertion_ids)
        assertions = civic.search_assertions_by_coordinates(query, search_mode='exact')
        assertion_ids = [x.id for x in assertions]
        assert set(assertion_ids) >= set(v600e_assertion_ids)

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

        query = CoordinateQuery('3', 10183706, 10183706, None, 'C')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('3', 10183706, 10183706, 'T', 'C')
        variants_single = civic.search_variants_by_coordinates(query, search_mode='exact')
        variants_bulk = civic.bulk_search_variants_by_coordinates([query], search_mode='exact')
        assert len(variants_single) == 1
        assert len(variants_bulk[query]) == 1
        assert hash(variants_single[0]) == variants_bulk[query][0].v_hash

        query = CoordinateQuery('3', 10183706, 10183706, '*', 'C')
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
        assert len(search_results[sorted_queries[0]]) == 19
        assert len(search_results[sorted_queries[1]]) >= 19

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
        assert len(search_results[sorted_queries[1]]) == 4

    def test_bulk_re_search_variants(self):
        sorted_queries = [
            CoordinateQuery('7', 140453136, 140453136),
            CoordinateQuery('7', 140453136, 140453137)
        ]
        search_results = civic.bulk_search_variants_by_coordinates(sorted_queries, search_mode='record_encompassing')
        assert len(search_results[sorted_queries[0]]) == 19
        assert len(search_results[sorted_queries[1]]) == 16

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
        assert len(search_results) == 1
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


class TestDrugs(object):

    def test_has_ncit_id(self, v600e_assertion):
        trametinib = v600e_assertion.drugs[0]
        assert trametinib.ncit_id == 'C77908'
        assert 'pubchem_id' not in trametinib.keys()

#warning logging tests
LOGGER = logging.getLogger(__name__)

def test_is_valid_for_vcf_warnings(caplog):
    fusion_variant = civic.get_variant_by_id(287)
    fusion_variant.is_valid_for_vcf(emit_warnings=True)
    assert "Variant 287 has a second set of coordinates. Skipping" in caplog.text

    incomplete_coordinates_variant = civic.get_variant_by_id(27)
    incomplete_coordinates_variant.is_valid_for_vcf(emit_warnings=True)
    assert "Incomplete coordinates for variant 27. Skipping." in caplog.text

    unsupported_var_bases_variant = civic.get_variant_by_id(613)
    unsupported_var_bases_variant.is_valid_for_vcf(emit_warnings=True)
    assert "Unsupported variant base(s) for variant 613. Skipping." in caplog.text

    #currently no case for unsupported ref bases
