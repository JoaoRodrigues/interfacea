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
Unit tests for FunctionalGroup classes and module.
"""

import copy
import os
import tempfile
import unittest

import interfacea as ia
import interfacea.functional_groups as fgs

TESTDIR = os.path.dirname(os.path.abspath(__file__))


class TestBasicFunctionalGroup(unittest.TestCase):
    """Simple test of an instance of the FunctionalGroup class.
    """

    def setUp(self):
        fpath = os.path.join(TESTDIR, 'data', 'mini.pdb')
        self.struct = ia.read(fpath)

    def test_assertions_and_exceptions(self):
        """Tests badly constructed FGs.
        """

        with self.assertRaisesRegexp(AssertionError, 'FunctionalGroup must have a name!'):
            _ = fgs.FunctionalGroup(charge=0,
                                    elements=[1, 1],
                                    bonds=[(0, 1)])

        with self.assertRaisesRegexp(AssertionError, 'FunctionalGroup must have a charge!'):
            _ = fgs.FunctionalGroup(name='test',
                                    elements=[1, 1],
                                    bonds=[(0, 1)])

        with self.assertRaisesRegexp(AssertionError, 'FunctionalGroup must have atoms!'):
            _ = fgs.FunctionalGroup(name='test',
                                    charge=0,
                                    bonds=[(0, 1)])

        with self.assertRaisesRegexp(AssertionError, 'FunctionalGroup must have bonds!'):
            _ = fgs.FunctionalGroup(name='test',
                                    charge=0,
                                    elements=[(1, 1)])

        with self.assertRaisesRegexp(AssertionError, 'Items in max_bonds must match elements'):
            _ = fgs.FunctionalGroup(name='test',
                                    charge=0,
                                    elements=[1, 1],
                                    bonds=[(0, 1)],
                                    max_bonds=[0])

        with self.assertRaises(fgs.FunctionalGroupError):
            _ = fgs.FunctionalGroup(name='isolated-element',
                                    charge=0,
                                    elements=[1, 1, 1],
                                    bonds=[(0, 1)])

        with self.assertRaises(fgs.FunctionalGroupError):
            _ = fgs.FunctionalGroup(name='extra-bond',
                                    charge=0,
                                    elements=[(1, 1)],
                                    bonds=[(0, 1), (0, 2)])

    def test_custom_FG(self):
        """Tests matching custom Functional Groups
        """

        # Neutral N Backbone (NH-CH-C-O)
        #
        #       H    H    O
        #       |    |    |
        #  X -- N -- C -- C --
        #            |
        #            X
        #
        elements = [7, 1, 6, 1, 6, 0, 8]
        bonds = [(0, 1), (0, 2), (2, 3), (2, 4), (2, 5), (4, 6)]
        max_bonds = [3, 1, 4, 1, 3, 4, 1]

        bb = fgs.FunctionalGroup(name='backbone',
                                 charge=0,
                                 elements=elements,
                                 bonds=bonds,
                                 max_bonds=max_bonds)

        matches = bb.search(self.struct)
        self.assertEqual(len(matches), 4)  # Same as no. of non-N-terminal residues

        # Any Backbone (NH-CH-C-O)
        #
        #       H    H    O
        #       |    |    |
        # 2X -- N -- C -- C --
        #            |
        #            X
        #
        elements = [7, 1, 6, 1, 6, 0, 8]
        bonds = [(0, 1), (0, 2), (2, 3), (2, 4), (2, 5), (4, 6)]

        bb = fgs.FunctionalGroup(name='any-backbone',
                                 charge=0,
                                 elements=elements,
                                 bonds=bonds)

        matches = bb.search(self.struct)
        self.assertEqual(len(matches), 6)  # Same as no. of residues

        # Terminal amine
        elements = [7, 1, 1, 1, 6]
        bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
        max_bonds = [4, 1, 1, 1, 4]
        nt = fgs.FunctionalGroup(name='terminal-amine',
                                 charge=1,
                                 elements=elements,
                                 bonds=bonds,
                                 max_bonds=max_bonds)

        matches = nt.search(self.struct)
        self.assertEqual(len(matches), 2)


class TestChargedFunctionalGroups(unittest.TestCase):
    """Test all charged functional groups.
    """

    def setUp(self):
        fpath = os.path.join(TESTDIR, 'data', 'mini.pdb')
        self.struct = ia.read(fpath)

    def test_findAnions(self):
        """Tests matching anionic functional groups
        """

        anions = fgs.anionic
        num_matches = 0
        for fg in anions:
            _fg = fg()
            m = _fg.search(self.struct)
            num_matches += len(m)

        self.assertEqual(num_matches, 2)  # Two aspartates

    def test_findCations(self):
        """Tests matching cationic functional groups
        """

        cations = fgs.cationic
        num_matches = 0
        for fg in cations:
            _fg = fg()
            m = _fg.search(self.struct)
            num_matches += len(m)

        self.assertEqual(num_matches, 4)  # Two N-terminal Arginines
