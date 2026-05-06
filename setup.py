from setuptools import setup, find_packages
from civicpy.__version__ import (
    __version__,
    __authors__,
    __author_email__,
    __description__,
    __url__,
)
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="civicpy",
    version=__version__,
    packages=find_packages(),
    url=__url__,
    license="MIT",
    author=", ".join(__authors__),
    author_email=__author_email__,
    description=__description__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    install_requires=[
        "requests",
        "obonet",
        "networkx",
        "pandas<=2.3.3",
        "Click",
        "vcfpy~=0.13.8",
        "pysam",
        "python-dateutil",
        "deprecation",
        "ga4gh.vrs @ git+https://github.com/ga4gh/vrs-python.git@d940af64126385b098365cc429a0313640a5d7f5",
        "ga4gh.cat_vrs @ git+https://github.com/ga4gh/cat-vrs-python.git@4c228e51a55b36752bbcc88012c3bb84a2299781",
        "ga4gh.va_spec @ git+https://github.com/ga4gh/va-spec-python.git@5ab6da0c09b6a9f01ddffa378f485a0d0aeeeb44",
    ],
    extras_require={
        "test": [
            "pytest==6.2.5",
            "pytest-cov==5.0.0",
            "attrs==22.1.0",
            "coveralls",
            "coverage<7.4.4",
            "deepdiff",
        ],
        "docs": ["sphinx", "sphinxjp.themes.basicstrap", "sphinxcontrib.programoutput"],
    },
    python_requires=">=3.10",
    entry_points={"console_scripts": ["civicpy=civicpy.cli:cli"]},
)
