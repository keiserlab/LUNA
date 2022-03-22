"""

LUNA: drug discovery toolkit.

If you have downloaded and uncompressed the LUNA source code, or fetched it from git,
for the simplest installation just type the command::

    python setup.py install

However, you would normally install the latest Biopython release from the PyPI archive
with::

    pip install luna

For more in-depth instructions, see the installation section of the
Biopython manual, linked to from:

http://luna-toolkit.readthedocs.io

Contributions, feature requests, and bug reports are greatly appreciated.
Please consult the issue tracker:

https://github.com/keiserlab/LUNA/issues

"""

from setuptools import setup, find_packages
from luna.version import version


__author__ = "Alexandre Fassio"
__date__ = "March 2022"
__maintainer__ = "Alexandre Fassio"
__email__ = "afassio@keiserlab.org"
__status__ = "Production"


def get_readme():
    with open('README.rst') as f:
        return f.read()


requirements = [
    'biopython==1.72',
    'colorlog',
    'matplotlib',
    'mmh3>=2.5.1',
    'networkx',
    'numpy',
    'pandas',
    'seaborn',
    'scipy',
    'xopen',
]

test_requirements = ["pytest", "mock"]

classifiers = ['Programming Language :: Python',
               'Programming Language :: Python :: 3.7',
               'Programming Language :: Python :: 3.8',
               'Programming Language :: Python :: 3.9',
               'Programming Language :: Python :: 3.10',
               'Programming Language :: Python :: 3.11',
               'License :: OSI Approved :: MIT License',
               'Operating System :: OS Independent',
               'Development Status :: 4 - Beta',
               'Intended Audience :: Science/Research',
               'Intended Audience :: Developers',
               'Intended Audience :: End Users/Desktop',
               'Topic :: Scientific/Engineering :: Bio-Informatics',
               'Topic :: Scientific/Engineering :: Chemistry',
               'Topic :: Software Development :: Libraries :: Python Modules'
               ]

setup(
    name='luna',
    version=version,
    packages=find_packages(exclude=['docs', 'tests']),
    python_requires='>=3.7, <4',
    install_requires=requirements,
    extras_require={
        'docs': ['sphinx >= 1.4', 'sphinx-rtd-theme']
    },
    include_package_data=True,
    description='LUNA: drug discovery toolkit',
    long_description=get_readme(),
    long_description_content_type='text/x-rst',
    url='https://github.com/keiserlab/LUNA',
    author='Alexandre Fassio',
    author_email='afassio@keiserlab.org',
    license='MIT',
    classifiers=classifiers,
    keywords='luna eifp fifp hifp ifp interaction fingerprint protein-ligand molecule docking',
    project_urls={
        'Bug Reports': 'https://github.com/keiserlab/LUNA/issues',
        'Source': 'https://github.com/keiserlab/LUNA/',
    },
    download_url='https://github.com/keiserlab/LUNA/tarball/v' + version,
    tests_require=test_requirements,

    entry_points={
        'console_scripts': [
            'run_luna = luna.run:main'
        ]
    }
)