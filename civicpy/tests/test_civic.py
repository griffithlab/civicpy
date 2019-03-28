import pytest
from civicpy import civic
from collections import defaultdict

ELEMENTS = [
    'Assertion'
]


def setup_module():
    try:
        civic.load_cache()
    except FileNotFoundError:
        pass
    civic.get_all_variants()


@pytest.fixture(scope="module", params=ELEMENTS)
def element(request):
    element_type = request.param
    return civic._get_elements_by_ids(element_type, [1])[0]


@pytest.fixture(scope="module")
def v600e():
    return civic.get_variant_by_id(12)


class TestGetFunctions(object):
    
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


class TestCoordinateSearch(object):

    def test_search_assertions(self):
        coordinates = {
            'chr': 7,
            'start': 140453136,
            'stop': 140453136,
            'alt': 'T'
        }
        assertions = civic.search_assertions_by_coordinates(coordinates)
        assertion_ids = [x.id for x in assertions]
        v600e_assertion_ids = (7, 10, 12, 20)
        v600k_assertion_ids = (11, 13)
        assert set(assertion_ids) == set(v600e_assertion_ids + v600k_assertion_ids)
        assertions = civic.search_assertions_by_coordinates(coordinates, search_mode='exact')
        assertion_ids = [x.id for x in assertions]
        assert set(assertion_ids) == set(v600e_assertion_ids)

    def test_bulk_search_variants(self):
        def coord_gen():
            coordinate_sets = [
                {
                    'chr': '7',
                    'start': 140453136,
                    'stop': 140453136,
                    'alt': 'T'
                },
                {
                    'chr': '7',
                    'start': 140453136,
                    'stop': 140453137,
                    'alt': 'TT'
                },
            ]
            for c in coordinate_sets:
                yield c
        gen = coord_gen()
        search_results = list(civic.bulk_search_variants_by_coordinates(gen))
        results_dict = defaultdict(list)
        for q, r in search_results:
            k = (q['chr'], q['start'], q['stop'], q['alt'])
            results_dict[k].append(civic.CACHE[r['v_hash']])
        active_variants = [v for v in results_dict[('7', 140453136, 140453136, 'T')] if len(v.evidence_items)]
        assert len(active_variants) >= 12
