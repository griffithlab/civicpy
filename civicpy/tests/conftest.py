from civicpy import LOCAL_CACHE_PATH
from civicpy.civic import load_cache


# Don't refresh stale cache data during tests
load_cache(LOCAL_CACHE_PATH, on_stale="ignore")
