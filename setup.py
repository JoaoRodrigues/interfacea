#!/usr/bin/env python

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
"""

#
# Borrows a lot from Biopython's setup.py script
#

import os
import sys

from setuptools import setup
# from setuptools import Command
from setuptools import Extension

from setuptools.command.install import install
from setuptools.command.build_py import build_py
from setuptools.command.build_ext import build_ext

__version__ = "0.1"

# Check for Python version
if sys.version_info[0] != 3:
    sys.stderr.write("interfacea requires Python 3.x. "
                     "Python %d.%d detected.\n" % sys.version_info[:2])
    sys.exit(1)

PACKAGES = [
    'interfacea',
    'interfacea.private',
]

EXTENSIONS = [
    Extension('interfacea.src.kdtrees',
              [os.path.join('interfacea', 'src', 'kdtrees.c')]),
]

REQUIRES = [
    'networkx',
    'numpy',
    'openmm',
    'pandas',
    'pdbfixer',
]

# For long description
with open("README.md", "rb") as handle:
    readme = handle.read().decode("ascii")

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


# class run_tests(Command):
#     """Run all of the tests for the package.
#     """

#     description = "Automatically run the test suite (in test/)"
#     user_options = []


setup(name='interfacea',
      version=__version__,
      author='Joao Rodrigues',
      author_email='j.p.g.l.m.rodrigues@gmail.com',
      url='https://github.com/joaorodrigues/interfacea',
      description='Open-source library to analyze the structure and energetics of biomolecular interfaces',
      long_description=readme,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: Freely Distributable',  # Change
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
          # "test": run_tests,
      },
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      install_requires=REQUIRES,
      )
