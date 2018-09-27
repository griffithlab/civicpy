import pytest
from civicpy import civic

ELEMENTS = [
    'Assertion'
]


@pytest.fixture(scope="module", params=ELEMENTS)
def element(request):
    element_type = request.param
    return civic._get_elements_by_ids(element_type, [1])[0]


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
            if isinstance(complex_value, civic.RecordAttribute):
                assert not complex_value._partial
