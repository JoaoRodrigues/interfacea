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
Unit tests for read() function
"""

import os
import unittest

import interfacea as ia

TESTDIR = os.path.dirname(os.path.abspath(__file__))


class TestRead(unittest.TestCase):

    def setUp(self):
        self.n_residues = 6
        self.n_atoms = 106
        self.n_chains = 2

    def test_missingFile(self):
        """Tests exception throwing when reading non-existent file.
        """

        fpath = os.path.join(TESTDIR, 'illusions', 'void.pdb')
        with self.assertRaises(ia.structure.StructureError):
            ia.read(fpath)

    def test_notSupportedExtension(self):
        """Tests exception throwing when reading file with unsupported extension.
        """

        fpath = os.path.join(TESTDIR, 'illusions', 'void.pdb')
        with self.assertRaises(ia.structure.StructureError):
            ia.read(fpath)

    def test_wrongFileType(self):
        """Tests exception throwing when reading file with wrong user-defined type.
        """

        fpath = os.path.join(TESTDIR, 'data', 'mini.pdb')
        with self.assertRaises(ia.structure.StructureError):
            ia.read(fpath, ftype='cif')

    def test_readPDB(self):
        """Tests reading/parsing a sample PDB file.
        """

        fpath = os.path.join(TESTDIR, 'data', 'mini.pdb')
        s = ia.read(fpath)
        top = s.topology

        self.assertEqual(top.getNumAtoms(), self.n_atoms)
        self.assertEqual(top.getNumResidues(), self.n_residues)
        self.assertEqual(top.getNumChains(), self.n_chains)

    def test_readCIF(self):
        """Tests reading/parsing a sample mmCIF file.
        """

        fpath = os.path.join(TESTDIR, 'data', 'mini.cif')
        s = ia.read(fpath)
        top = s.topology

        self.assertEqual(top.getNumAtoms(), self.n_atoms)
        self.assertEqual(top.getNumResidues(), self.n_residues)
        self.assertEqual(top.getNumChains(), self.n_chains)
