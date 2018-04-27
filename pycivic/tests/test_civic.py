import pytest
from pycivic import civic

ELEMENTS = [
    'Assertion'
]

UNMARKED_PLURALS = set(
    'Evidence'
)


@pytest.fixture(scope="module", params=ELEMENTS)
def element(request):
    element_type = request.param
    t = element_type.lower()
    if element_type in UNMARKED_PLURALS:
        f = getattr(civic.MODULE, f'get_{t}')
    else:
        f = getattr(civic.MODULE, f'get_{t}s')
    return f([1])[0]


class TestGetFunctions(object):
    
    def test_get_assertions(self):
        test_ids = [1, 2, 3]
        results = civic.get_assertions(test_ids)
        assert len(results) == 3


class TestCivicRecord(object):

    def test_module(self):
        assert str(type(civic.MODULE)) == "<class 'module'>"


class TestElements(object):

    def test_attribute_fail(self, element):
        with pytest.raises(AttributeError):
            element.foo

    def test_completeness(self, element):
        for complex_field in element.COMPLEX_FIELDS:
            complex_value = getattr(element, complex_field)
            if not complex_value:
                continue
            if isinstance(complex_value, list):
                complex_value = complex_value[0]
            assert complex_value.partial == False
            if complex_value.type in ['variant', 'evidence', 'gene']:
                assert complex_value.description
