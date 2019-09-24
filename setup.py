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
        'pytest',
        'pytest-cov',
        'python-coveralls',
        'requests',
        'obonet',
        'networkx',
        'pandas',
        'Click'
    ],
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'civicpy=civicpy.cli:cli'
        ]
    },
)
