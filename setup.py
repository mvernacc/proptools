'''A setuptools based installer for proptools.

Based on https://github.com/pypa/sampleproject/blob/master/setup.py

Matt Vernacchia
proptools
2016 Sept 21
'''

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

INSTALL_REQUIRES = [
    'numpy',
    'matplotlib',
    'scikit-aero',
    'scipy',
    ]
TEST_REQUIRES = [
        'pytest',
        'coverage',
        'pytest-cov',
    ]
DOCS_REQUIRES = [
    'sphinx_rtd_theme',
    ]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='proptools-rocket',

    version='0.0.2',

    description='Rocket propulsion design calculation tools.',

    author='Matt Vernacchia',
    author_email='mvernacc@mit.edu',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.7',
    ],

    # What does your project relate to?
    keywords='rocket propulsion engineering aerospace',

    long_description=long_description,
    long_description_content_type="text/markdown",

    project_urls={
        'Documentation': 'https://proptools.readthedocs.io/en/latest/',
        'Source Code': 'https://github.com/mvernacc/proptools',
    },

    install_requires=INSTALL_REQUIRES,
    extras_require={
        'test': TEST_REQUIRES + INSTALL_REQUIRES,
        'docs': DOCS_REQUIRES + INSTALL_REQUIRES,
    },

    packages=find_packages(),

    scripts=[],
)
