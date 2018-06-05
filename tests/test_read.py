#!/usr/bin/env python

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
