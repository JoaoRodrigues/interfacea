#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2020 João Pedro Rodrigues
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
Unit tests for io.pdb module classes.
"""

import pathlib

import interfacea
from interfacea.io import pdb


class TestPDBReader:
    """Tests for PDBReader class."""

    def setup(self):
        """Pre-test setup."""

        rootdir = pathlib.Path('.')  # rootdir
        self.datadir = rootdir / 'tests' / 'data'

    def test_read_pdb(self):
        """Load and parse a PDB file."""

        data = pdb.PDBReader(self.datadir / 'mini.pdb')
