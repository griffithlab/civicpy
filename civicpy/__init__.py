from .__version__ import __version__
from pathlib import Path
import os

REMOTE_CACHE_URL = os.getenv('CIVICPY_REMOTE_CACHE_URL', False) or \
    "https://civicdb.org/downloads/nightly/nightly-civicpy_cache.pkl"
LOCAL_CACHE_PATH = os.getenv('CIVICPY_CACHE_FILE', False) or \
    str(Path.home() / '.civicpy' / 'cache.pkl')
CACHE_TIMEOUT_DAYS = os.getenv('CIVICPY_CACHE_TIMEOUT_DAYS', False) or 7
# We convert the CACHE_TIMEOUT_DAYS to int
# if not possible we use the default value
try:
    CACHE_TIMEOUT_DAYS = int(CACHE_TIMEOUT_DAYS)
except TypeError:
    CACHE_TIMEOUT_DAYS = 7

TEST_CACHE_PATH = str(Path(__file__).parent / 'data' / 'test_cache.pkl')


def version():
    return __version__
