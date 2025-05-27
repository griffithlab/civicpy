from setuptools import setup, find_packages
from civicpy.__version__ import __version__, __authors__, __author_email__, __description__, __url__
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='civicpy',
    version=__version__,
    packages=find_packages(),
    url=__url__,
    license='MIT',
    author=', '.join(__authors__),
    author_email=__author_email__,
    description=__description__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
    ],
    install_requires=[
        'requests',
        'obonet',
        'networkx',
        'pandas',
        'Click',
        'vcfpy',
        'pysam',
        'backports-datetime-fromisoformat',
        'deprecation',
        'ga4gh.vrs',
        'ga4gh.cat_vrs',
        'ga4gh.va_spec~=0.3.0'
    ],
    extras_require={
        'test': [
            'pytest==6.2.5',
            'pytest-cov==5.0.0',
            'attrs==22.1.0',
            'coveralls',
            'coverage<7.4.4',
            'deepdiff'
        ],
        'docs': [
            'sphinx',
            'sphinxjp.themes.basicstrap',
            'sphinxcontrib.programoutput'
        ]
    },
    python_requires='>=3.10',
    entry_points={
        'console_scripts': [
            'civicpy=civicpy.cli:cli'
        ]
    },
)
