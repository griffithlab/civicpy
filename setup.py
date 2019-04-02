from setuptools import setup
from civicpy.__version__ import __version__, __author__, __author_email__, __description__, __url__

setup(
    name='civicpy',
    version=__version__,
    packages=['civicpy'],
    url=__url__,
    license='MIT',
    author=__author__,
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
        'requests',
        'obonet',
        'networkx',
        'pandas'
    ],
    python_requires='>=3.7',
    entry_points={},
)
