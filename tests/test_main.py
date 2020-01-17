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
Unit tests for main module functions (__init__.py)
"""

import logging
import os
import unittest

import interfacea as ia
from interfacea.core.Structure import StructureError

TESTDIR = os.path.dirname(os.path.abspath(__file__))


class TestLogger(unittest.TestCase):

    def tearDown(self):
        """Reset log levels"""
        ia.set_log_level(level='silent')

    def test_default_log_level(self):
        """Tests the default setup of the logging framework"""

        root_logger = logging.getLogger()
        log_level = root_logger.getEffectiveLevel()
        self.assertEqual(logging.CRITICAL, log_level)

    def test_change_log_level(self):
        """Tests changing the log level to verbose"""

        ia.set_log_level(level='verbose')

        root_logger = logging.getLogger()
        log_level = root_logger.getEffectiveLevel()
        self.assertEqual(logging.INFO, log_level)

    def test_change_log_error(self):
        """Tests changing the log level to an unsupported value"""

        with self.assertRaises(ValueError):
            ia.set_log_level(level='bananas')


class TestRandomSeed(unittest.TestCase):

    def test_check_random_seed(self):
        """Tests access to the global random seed variable"""

        try:
            _ = ia.RANDOM_SEED
        except Exception as e:
            self.fail(f'Exception raised when accessing RANDOM_SEED {e}')

    def test_change_random_seed(self):
        """Tests changing the random seed variable"""

        new_seed = 2907
        ia.set_random_seed(new_seed)
        self.assertEqual(new_seed, ia.RANDOM_SEED)

    def test_invalid_random_seed_str(self):
        """Tests changing the random seed variable to a string"""

        new_seed = 'abc'
        with self.assertRaises(TypeError):
            ia.set_random_seed(new_seed)

    def test_invalid_random_seed_neg(self):
        """Tests changing the random seed variable to a negative number"""

        new_seed = -10
        with self.assertRaises(TypeError):
            ia.set_random_seed(new_seed)


class TestRead(unittest.TestCase):

    def setUp(self):

        self.mini_residues = 6
        self.mini_atoms = 106
        self.mini_chains = 2

        self.noexists = os.path.join(TESTDIR, 'data', 'void.pdb')  # missing
        self.fake_ext = os.path.join(TESTDIR, 'data', 'mini.abc')
        self.mini_pdb = os.path.join(TESTDIR, 'data', 'mini.pdb')
        self.mini_cif = os.path.join(TESTDIR, 'data', 'mini.cif')

    def test_missingFile(self):
        """Exception thrown when reading non-existent file."""

        with self.assertRaises(StructureError):
            ia.read(self.noexists)

    def test_notSupportedExtension(self):
        """Exception thrown when reading file with unsupported extension."""

        with self.assertRaises(StructureError):
            ia.read(self.fake_ext)

    def test_read_pdb_as_cif(self):
        """Test exception thrown when reading mmCIF file as PDB."""

        with self.assertRaises(StructureError):
            ia.read(self.mini_pdb, ftype='cif')

    def test_read_cif_as_pdb(self):
        """Test exception thrown when reading PDB file as mmCIF."""

        with self.assertRaises(StructureError):
            ia.read(self.mini_cif, ftype='pdb')

    def test_read_pdb(self):
        """Tests reading/parsing a sample PDB file.
        """

        s = ia.read(self.mini_pdb)
        top = s.topology

        self.assertEqual(top.getNumAtoms(), self.mini_atoms)
        self.assertEqual(top.getNumResidues(), self.mini_residues)
        self.assertEqual(top.getNumChains(), self.mini_chains)

    def test_read_cif(self):
        """Tests reading/parsing a sample mmCIF file.
        """

        s = ia.read(self.mini_cif)
        top = s.topology

        self.assertEqual(top.getNumAtoms(), self.mini_atoms)
        self.assertEqual(top.getNumResidues(), self.mini_residues)
        self.assertEqual(top.getNumChains(), self.mini_chains)


if __name__ == "__main__":

    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
