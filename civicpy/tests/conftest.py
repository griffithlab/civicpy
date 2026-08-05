from civicpy.civic import load_cache, get_molecular_profile_by_id
import pytest
from ga4gh.vrs.models import Allele, CopyNumberChange


@pytest.fixture(scope="module")
def braf_v600e_vrs():
    params = {
        "name": "BRAF V600E",
        "id": "ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
        "type": "Allele",
        "digest": "j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
        "location": {
            "id": "ga4gh:SL.t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
            "type": "SequenceLocation",
            "digest": "t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
            },
            "start": 599,
            "end": 600,
            "sequence": "V",
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
    }
    return Allele(**params)


@pytest.fixture(scope="session")
def braf_amplification_vrs():
    params = {
        "name": "BRAF Amplification",
        "id": "ga4gh:CX.h7unj-f_djER28-h2Q6Prvo3C90O4d3M",
        "type": "CopyNumberChange",
        "digest": "h7unj-f_djER28-h2Q6Prvo3C90O4d3M",
        "location": {
            "id": "ga4gh:SL.0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
            "type": "SequenceLocation",
            "digest": "0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            },
            "start": 140713327,
            "end": 140924929,
        },
        "copyChange": "high-level gain",
    }
    return CopyNumberChange(**params)


@pytest.fixture(scope="module")
def v600e_mp():
    return get_molecular_profile_by_id(12)
