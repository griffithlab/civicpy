import pytest
from pycivic import civic


class TestGetFunctions(object):
    
    def test_get_assertions(self):
        test_ids = [1, 2, 3]
        results = civic.get_assertions(test_ids)
        assert len(results) == 3
        assertion_1 = [x for x in results if x.id == 1][0]
        assert assertion_1.variant.id == 306


class TestCivicRecord(object):

    def test_module(self):
        assert str(type(civic.MODULE)) == "<class 'module'>"
