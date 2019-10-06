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
Unit tests for Structure class and methods
"""

import copy
import os
import tempfile
import unittest
import warnings

import interfacea as ia

TESTDIR = os.path.dirname(os.path.abspath(__file__))

warnings.simplefilter('ignore')  # blanket warning ignore filter


class TestReader(unittest.TestCase):
    """Tests for read() method and associated functions.
    """

    def setUp(self):
        fpath = os.path.join(TESTDIR, 'data', 'protein_ligand.pdb')
        self.struct = ia.read(fpath)

    def test_read_ligand(self):
        """Test writing a Structure with a ligand (not in topology).
        """

        struct = self.struct
        # Check DRG residue
        drg = [r for r in struct.topology.residues() if r.name == 'DRG'][0]
        n_bonds = len(list(drg.bonds()))
        self.assertEqual(n_bonds, 26)

    def test_make_residue_graph(self):
        """Tests creating a nx.Graph from residue topology.
        """

        struct = self.struct
        reslist = [r for r in struct.topology.residues()]
        for residue in reslist:
            g = residue._g

            n_atoms = len(list(residue.atoms()))
            n_nodes = g.number_of_nodes()
            self.assertEqual(n_atoms, n_nodes)

            n_bonds = len(list(residue.internal_bonds()))
            n_edges = g.number_of_edges()
            self.assertEqual(n_bonds, n_edges)

            # Check C-alpha has 4 edges
            if residue.name == 'DRG':
                continue

            idx = [idx for idx, a in enumerate(residue.atoms())
                   if a.name == 'CA'][0]
            n_node_edges = len(g[idx])
            self.assertEqual(n_node_edges, 4)




class TestTermini(unittest.TestCase):
    """Tests for add_termini() method
    """

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


class TestAddAtoms(unittest.TestCase):
    """Tests for add_missing_atoms() method
    """

    def setUp(self):
        self.s_allatom = ia.read(os.path.join(TESTDIR, 'data', 'mini.pdb'))
        self.s_missing = ia.read(os.path.join(TESTDIR, 'data', 'mini_incomplete.pdb'))

    def test_addToCompleteStructure(self):
        """Test adding atoms to structure not missing any.
        """

        s = copy.deepcopy(self.s_allatom)  # copy initial structure
        s.add_missing_atoms()

        n_ini_atoms = self.s_allatom.topology.getNumAtoms()
        n_end_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_ini_atoms, n_end_atoms)

    def test_addToIncompleteStructure(self):
        """Test adding atoms to structure missing 1 side-chain atom.
        """

        s = copy.deepcopy(self.s_missing)  # copy initial structure
        s.add_missing_atoms()

        n_ini_atoms = self.s_missing.topology.getNumAtoms()
        n_end_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_end_atoms - n_ini_atoms, 1)


class TestMutate(unittest.TestCase):
    """Tests for mutate() method
    """

    def setUp(self):
        self.struct = ia.read(os.path.join(TESTDIR, 'data', 'mini_noH.pdb'))

    def test_single_mutation(self):
        """Test performing a single mutation.
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        s.mutate([('A-ASN-1', 'ALA')])

        self.assertEqual(list(s.topology.residues())[0].name, 'ALA')

        n_ini_atoms = self.struct.topology.getNumAtoms()
        n_end_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_end_atoms - n_ini_atoms, -3)

    def test_multiple_mutations(self):
        """Test performing multiple mutations at once.
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        s.mutate([('A-ASN-1', 'ALA'), ('B-ALA-6', 'PHE')])

        reslist = [r.name for r in s.topology.residues()]
        self.assertEqual(reslist[0], 'ALA')
        self.assertEqual(reslist[-1], 'PHE')

        n_ini_atoms = self.struct.topology.getNumAtoms()
        n_end_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_end_atoms - n_ini_atoms, 3)

    def test_mutateNonExistentChain(self):
        """Test throwing exception on mutating a residue on a non-existent chain.
        """

        s = copy.deepcopy(self.struct)  # copy initial structure
        with self.assertRaises(ia.structure.StructureError):
            s.mutate([('C-ASN-1', 'ALA')])


class TestWriter(unittest.TestCase):
    """Tests for write() method.
    """

    def setUp(self):
        self.pdb = os.path.join(TESTDIR, 'data', 'mini.pdb')
        self.cif = os.path.join(TESTDIR, 'data', 'mini.cif')

        self.struct = ia.read(self.pdb)

    def test_write_PDB(self):
        """Test writing a Structure in PDB format.
        """

        s = self.struct

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        s.write(filename, ftype='pdb', overwrite=True)
        with open(filename, 'r') as handle:
            new = handle.readlines()

        os.remove(filename)

        with open(self.pdb) as handle:
            ori = handle.readlines()

        self.assertEqual(new[1:], ori[1:])  # remove remark line w/ date.

# Broken in latest version of openMM
#     def test_write_mmCIF(self):
#         """Test writing a Structure in mmCIF format.
#         """

#         s = self.struct

#         filenumber, filename = tempfile.mkstemp()
#         os.close(filenumber)

#         s.write(filename, ftype='cif', overwrite=True)
#         with open(filename, 'r') as handle:
#             new = handle.readlines()

#         # os.remove(filename)

#         with open(self.cif, 'r') as handle:
#             ori = handle.readlines()

#         with open('test.cif', 'w') as handle:
#             handle.write(''.join(new))

#         self.assertEqual(new[2:], ori[2:])  # remove lines w/ date.

    def test_writeToExistingFile(self):
        """Test writing a Structure to an existing file.
        """

        s = self.struct
        with self.assertRaises(OSError):
            with tempfile.NamedTemporaryFile() as handle:
                s.write(handle.name, ftype='pdb')
