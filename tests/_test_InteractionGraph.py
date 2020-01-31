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
Unit tests for Structure class and methods
"""

import os
import unittest
import warnings

from simtk.openmm.app import Residue as _Residue, Chain as _Chain

import interfacea.core.InteractionGraph as ig

TESTDIR = os.path.dirname(os.path.abspath(__file__))

warnings.simplefilter('ignore')  # blanket warning ignore filter


class TestInteractionGraphSimple(unittest.TestCase):
    """Tests for InteractionGraph class - Simple Cases.
    """

    def setUp(self):

        # Dummy Residue Mock
        c = _Chain(0, None, 'A')
        self.r1 = _Residue('ARG', 0, c, 1, ' ')
        self.r2 = _Residue('GLU', 0, c, 2, ' ')

        # Dummy Graphs
        self.g1 = ig.InteractionGraph(name='genericGraph')
        self.g1.add_residue(self.r1)

        self.g2 = ig.InteractionGraph(name='genericGraph')
        self.g2.add_residue(self.r1)
        self.g2.add_residue(self.r2)

        self.g2e = ig.InteractionGraph(name='genericGraph')
        self.g2e.add_residue(self.r1)
        self.g2e.add_residue(self.r2)
        self.g2e.add_interaction(self.r1, self.r2, {'type': 'geometric'})

    # __contains__
    def test_contains_residue(self):
        """Tests __contains__ method for single Residues"""
        self.assertTrue(self.r1 in self.g1)
        self.assertFalse(self.r2 in self.g1)

    def test_contains_interaction(self):
        """Tests __contains__ method for interactions"""
        self.assertTrue((self.r1, self.r2) in self.g2e)
        self.assertTrue((self.r2, self.r1) in self.g2e)  # invariant to order
        self.assertFalse((self.r1, self.r2) in self.g2)

    def test_contains_interaction_with_attributes(self):
        """Tests __contains__ method for interactions with attributes"""
        self.assertTrue((self.r2, self.r1, {'type': 'geometric'}) in self.g2e)
        self.assertFalse((self.r2, self.r1, {'type': 'energy'}) in self.g2e)

    # # Node Tests
    def test_add_new_residue(self):
        """Test adding a mock residue to the InteractionGraph.
        """

        g = ig.InteractionGraph(name='genericGraph')
        idx = g.add_residue(self.r1)
        self.assertEqual(idx, 1)

    def test_add_existing_residue(self):
        """Test adding an existing residue to InteractionGraph.
        """

        idx = self.g1.add_residue(self.r1)
        self.assertEqual(idx, 1)

    def test_add_existing_residue_strict(self):
        """Test exception when adding an existing residue to the InteractionGraph.
        """

        with self.assertRaises(ig.DuplicatedNodeException):
            self.g1.add_residue(self.r1, ignore_existing=False)

    def test_add_another_residue(self):
        """Test adding a new residue to the InteractionGraph.
        """

        idx = self.g1.add_residue(self.r2)
        self.assertEqual(idx, 2)

    # Edge Tests
    def test_add_new_edge_simple(self):
        """Test adding new interaction between existing residues.
        """

        idx = self.g2.add_interaction(self.r1, self.r2)
        self.assertEqual(idx, 0)  # edges start at 0

        # Make sure edge was properly added
        self.assertTrue((self.r1, self.r2) in self.g2)

    def test_add_new_edge_kw(self):
        """Test adding interaction with attributes between existing residues.
        """

        idx = self.g2.add_interaction(self.r1, self.r2, {'type': 'geometric'})
        self.assertEqual(idx, 0)

        # Retrieve edge data and compare
        expected = {'type': 'geometric'}
        real = self.g2._g.get_edge_data(1, 2)[idx]
        self.assertEqual(real, expected)

    def test_add_two_edges(self):
        """Test adding two interactions between existing residues.
        """

        idx = self.g2.add_interaction(self.r1, self.r2, {})
        self.assertEqual(idx, 0)
        idx = self.g2.add_interaction(self.r1, self.r2, {'type':'geometric'})
        self.assertEqual(idx, 1)


if __name__ == "__main__":

    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
