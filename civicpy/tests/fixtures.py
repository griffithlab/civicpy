import pytest
from civicpy import civic


@pytest.fixture()
def v600e():
    return civic.get_variant_by_id(12)