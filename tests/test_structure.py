#!/usr/bin/env python

"""
Unit tests for Structure class and methods
"""

import copy
import os
import unittest

import interfacea as ia

TESTDIR = os.path.dirname(os.path.abspath(__file__))


class TestStructure(unittest.TestCase):

    def setUp(self):
        fpath = os.path.join(TESTDIR, 'data', 'mini.pdb')
        self.struct = ia.read(fpath)

    def test_wrongTerminiType(self):
        """Test exception throwing when user defines unknown cap.
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        with self.assertRaises(ia.structure.StructureError):
            s.add_termini(ends=[('ACE', 'NME'), ('ABC', 'NME')])

    def test_addDefaultTermini(self):
        """Test default termini caps (ACE/NME).
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        s.add_termini()

        for chain in s.topology.chains():
            reslist = [r.name for r in chain.residues()]
            self.assertEqual(reslist[0], 'ACE')
            self.assertEqual(reslist[-1], 'NME')

        # Assert method did not add any other atoms
        n_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_atoms, 116)

    def test_wrongNumEnds(self):
        """Test exception throwing when caps do not match number of chains.
        """

        s = copy.deepcopy(self.struct)  # copy initial structure

        with self.assertRaises(ia.structure.StructureError):
            s.add_termini(ends=[(None, None)])

        with self.assertRaises(ia.structure.StructureError):
            s.add_termini(ends=[(None, None), (None, None), (None, None)])

    def test_addUncappedTermini(self):
        """Test uncapped termini (add OXT atoms).
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        s.add_termini(ends=[(None, None), (None, None)])

        for chain in s.topology.chains():
            terminal_res = [r for r in chain.residues()][-1]
            atoms = {a.name for a in terminal_res.atoms()}
            self.assertIn('OXT', atoms)

        # Assert method did not add any other atoms
        n_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_atoms, 108)

    def test_addMixedTermini(self):
        """Test mixed (uncapped/capped) termini.
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        s.add_termini(ends=[('ACE', None), (None, 'NME')])

        for chain in s.topology.chains():
            if chain.index == 0:
                reslist = [r for r in chain.residues()]
                c_ter_atoms = {a.name for a in reslist[-1].atoms()}
                self.assertEqual(reslist[0].name, 'ACE')
                self.assertNotEqual(reslist[-1].name, 'NME')
                self.assertIn('OXT', c_ter_atoms)

            elif chain.index == 1:
                reslist = [r for r in chain.residues()]
                self.assertNotEqual(reslist[0].name, 'ACE')
                self.assertEqual(reslist[-1].name, 'NME')

        # Assert method did not add any other atoms
        n_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_atoms, 112)