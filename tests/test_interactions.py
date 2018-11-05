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
Unit tests for interactions.py classes and module.
"""

import os
import unittest

import interfacea as ia

TESTDIR = os.path.dirname(os.path.abspath(__file__))


class TestInteractionAnalyzer(unittest.TestCase):
    """Test suite for the InteractionAnalyzer class.
    """

    def setUp(self):
        pp_path = os.path.join(TESTDIR, 'data', 'protein_protein.cif')
        pl_path = os.path.join(TESTDIR, 'data', 'protein_ligand.pdb')

        pp_complex = ia.read(pp_path)
        pl_complex = ia.read(pl_path)

        self.pp_ia = ia.InteractionAnalyzer(pp_complex)
        self.pl_ia = ia.InteractionAnalyzer(pl_complex)

    # Test Finders
    def test_find_cations(self):
        """testing find_cations()
        """

        pp_ia = self.pp_ia
        pp_ia.cations = None  # clear just in case
        pp_ia.find_cations()

        n_res_cations = len(pp_ia.cations)  # should be 24
        n_cations = sum(map(len, pp_ia.cations.values()))  # should be 25
        self.assertEqual(n_res_cations, 24)
        self.assertEqual(n_cations, 25)

    def test_find_anions(self):
        """testing find_anions()
        """

        pp_ia = self.pp_ia
        pp_ia.anions = None
        pp_ia.find_anions()

        n_res_anions = len(pp_ia.anions)
        n_anions = sum(map(len, pp_ia.anions.values()))
        self.assertEqual(n_res_anions, 28)
        self.assertEqual(n_anions, 28)

    def test_find_hydrophobic(self):
        """testing find_hydrophobics()
        """

        pp_ia = self.pp_ia
        pp_ia.hydrophobics = None
        pp_ia.find_hydrophobics()

        n_res = len(pp_ia.hydrophobics)
        n_groups = sum(map(len, pp_ia.hydrophobics.values()))
        self.assertEqual(n_res, 166)
        self.assertEqual(n_groups, 166)
