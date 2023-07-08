#!/usr/bin/env python

from setuptools import setup, find_packages
import os
from os import path

here = path.abspath(path.dirname(__file__))

with open('README.md', 'r') as rmf:
    readme = rmf.read()

############

setup(
    name='cx_analysis',
    version='0.0.1',
    author='Nicholas Chua',
    author_email='nchua@flatironinstitute.org',
    url='https://github.com/nicholasjchua/cx-analysis',
    description='Analysis scripts for connectomics research.',
    long_description=readme,
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
        'Programming Language :: Python :: 3',
    ],
    keywords='neuroscience analysis connectomics catmaid',
    packages=find_packages(),
    install_requires=['']
)
