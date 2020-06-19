from setuptools import setup
from civicpy.__version__ import __version__, __authors__, __author_email__, __description__, __url__

setup(
    name='civicpy',
    version=__version__,
    packages=['civicpy'],
    url=__url__,
    license='MIT',
    author=', '.join(__authors__),
    author_email=__author_email__,
    description=__description__,
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
    ],
    extras_require={
        'test': [
            'pytest==4.1.0',
            'pytest-cov==2.9.0',
            'attrs==18.2.0',
            'python-coveralls',
            'coverage<5.0',
        ],
        'docs': [
            'sphinx',
            'sphinxjp.themes.basicstrap',
            'sphinxcontrib.programoutput'
        ]
    },
    python_requires='>=3.5',
    entry_points={
        'console_scripts': [
            'civicpy=civicpy.cli:cli'
        ]
    },
)
