from setuptools import setup

setup(
    name='civicpy',
    version='0.0.1.a4',
    packages=['civicpy'],
    url='https://github.com/griffithlab/civicpy',
    license='MIT',
    author='Alex H. Wagner',
    author_email='awagner24@wustl.edu',
    description='a python wrapper and tools for the civic API',
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
        'requests'
    ],
    python_requires='~=3.6',
    entry_points={},

)
