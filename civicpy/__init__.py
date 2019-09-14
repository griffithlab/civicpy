from .__version__ import __version__
from pathlib import Path
import os

REMOTE_CACHE_URL = os.getenv('CIVICPY_REMOTE_CACHE_URL', False) or \
                   "https://civicdb.org/downloads/civicpy_cache.pkl"
LOCAL_CACHE_PATH = os.getenv('CIVICPY_CACHE_FILE', False) or \
                   str(Path.home() / '.civicpy' / 'cache.pkl')

def version():
    return __version__
