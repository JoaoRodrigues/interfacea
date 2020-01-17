#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018 João Pedro Rodrigues
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Setuptools-based setup/install script for interfacea

This uses setuptools, the standard python mechanism for installing packages.
After downloading and uncompressing the source code, or cloning it from git,
type the following command:

    python setup.py install

As a developer, to test your changes without installing the package please use:
    python setup.py develop
"""

#
# Borrows a lot from Biopython's setup.py script
#

import os
from pathlib import Path
import sys

from setuptools import setup
from setuptools import Extension

from setuptools.command.install import install
from setuptools.command.build_py import build_py
from setuptools.command.build_ext import build_ext

# Check for Python version
if sys.version_info < (3, 6):
    sys.stderr.write("interfacea requires Python 3.6+. "
                     "Python %d.%d detected.\n" % sys.version_info[:2])
    sys.exit(1)


# Get version from code
def get_version():
    """Read the variable version from interfacea/_version.py"""

    f = Path('./interfacea') / '_version.py'
    contents = f.read_text()
    version = contents.strip().split()[-1]
    return version[1:-1]


# Get long description
def get_long_description():
    """Reads the contents of the README file"""
    contents = Path('README.rst').read_text()
    return contents


PACKAGES = [
    'interfacea',
    'interfacea.core',
    'interfacea.analyzers',
    'interfacea.chemistry',
]

EXTENSIONS = [
    Extension('interfacea.src.kdtree.kdtrees',
              [os.path.join('interfacea', 'src', 'kdtree', 'kdtrees.c')]),
]

REQUIRES = [
    'networkx',
    'numpy',
    'openmm',
    'pdbfixer',
]

#
# Install/Build/Test classes
#


class install_library(install):
    def run(self):
        """Run the installation."""
        install.run(self)


class build_py_modules(build_py):
    def run(self):
        """Run the build."""
        build_py.run(self)


class build_extensions(build_ext):
    def run(self):
        """Run the build."""
        build_ext.run(self)


setup(name='interfacea',
      version=get_version(),
      author='Joao Rodrigues',
      author_email='j.p.g.l.m.rodrigues@gmail.com',
      url='https://github.com/joaorodrigues/interfacea',
    #   download_url=''  # TODO: add url to release using version string.
      description='Open-source library to analyze the structure and energetics of protein interfaces',
      long_description=get_long_description(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: Apache Software License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Bioinformatics',
          'Topic :: Scientific/Engineering :: Computational Biology',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      cmdclass={
          "install": install_library,
          "build_py": build_py_modules,
          "build_ext": build_extensions,
      },
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      install_requires=REQUIRES,
      )
