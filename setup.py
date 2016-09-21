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

setup(
    name='proptools',
    
    version='0.0.0a0',

    description='Rocket propulsion design calculation tools.',

    author='Matt Vernacchia',
    author_email='mvernacc@mit.edu.',

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
        'Programming Language :: Python :: 2.7',
    ],

     # What does your project relate to?
    keywords='rocket propulsion engineering aerospace',


    install_requires=['numpy'],

    packages=find_packages(),

    scripts=[],

)
