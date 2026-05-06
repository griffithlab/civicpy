from .__version__ import __version__
from pathlib import Path
import os
import requests
import subprocess
from packaging.version import Version
import logging
from datetime import datetime


#get the currently installed civicpy version and warn if it is out of date
result = subprocess.run(["pip", "index", "versions", "civicpy"], capture_output=True, text=True)
latest_version = result.stdout.split("LATEST:")[1].strip()
if Version(__version__) < Version(latest_version):
    logging.warning("The installed civicpy version is out of date.")
    logging.warning(f"Latest version: {latest_version}")
    logging.warning(f"Installed version: {__version__}")

REMOTE_CACHE_URL = os.getenv('CIVICPY_REMOTE_CACHE_URL', False)
fallback_cache_timeout_days = 7

if REMOTE_CACHE_URL:
    logging.warning(f"Using REMOTE_CACHE_URL: {REMOTE_CACHE_URL}")
else:
    versioned_url = f"https://civicdb.org/downloads/nightly/nightly-civicpy_cache-{__version__}.pkl"
    generic_url = "https://civicdb.org/downloads/nightly/nightly-civicpy_cache.pkl"
    use_versioned_url = False
    response = requests.head(versioned_url, timeout=5)
    # Checks if the status code is in the 200-399 range
    if response.status_code == requests.codes.ok:
        use_versioned_url = True
        last_modified = response.headers.get('Last-Modified')
    else:
        logging.warning(f"Remote cache URL matching your installed civicpy version not found.")

    if use_versioned_url:
        REMOTE_CACHE_URL = versioned_url
        logging.warning(f"Using remote cache URL matching your installed civicpy version: {versioned_url}")
        delta = (datetime.today()-datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S GMT")).days
        logging.warning(f"This remote cache was created on {last_modified} and is {delta} days old.")
        if Version(__version__) < Version(latest_version):
            logging.warning("If a newer cache is desired, please update CIViCpy to the latest version and try again.")
        if delta >=7:
            fallback_cache_timeout_days = delta + 1
    else:
        REMOTE_CACHE_URL = generic_url
        logging.warning(f"Using generic remote cache URL: {generic_url}")
        if Version(__version__) < Version(latest_version):
            logging.warning(f"This cache might not be compatible with your civicpy version. If loading fails, please update CIViCpy to the latest version and try again.")

LOCAL_CACHE_PATH = os.getenv('CIVICPY_CACHE_FILE', False) or \
    str(Path.home() / '.civicpy' / 'cache.pkl')

CACHE_TIMEOUT_DAYS = os.getenv('CIVICPY_CACHE_TIMEOUT_DAYS', False) or fallback_cache_timeout_days
# We convert the CACHE_TIMEOUT_DAYS to int
# if not possible we use the default value
try:
    CACHE_TIMEOUT_DAYS = int(CACHE_TIMEOUT_DAYS)
except TypeError:
    CACHE_TIMEOUT_DAYS = fallback_cache_timeout_days
logging.warning(f"Setting cache timeout days to {CACHE_TIMEOUT_DAYS}.")
