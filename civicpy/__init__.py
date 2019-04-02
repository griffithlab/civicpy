from .__version__ import __version__
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
DATA_ROOT = PROJECT_ROOT / 'data'

if not DATA_ROOT.exists():
    DATA_ROOT.mkdir()


def version():
    return __version__