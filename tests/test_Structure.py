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

import pathlib
import shutil
import tempfile
import unittest
import warnings

import simtk.openmm.app as app

import interfacea as ia
from interfacea.core.Structure import Structure
from interfacea.core.Structure import (
    StructureError,
    StructureWriteError
)


TESTDIR = pathlib.Path(__file__).resolve().parent


class TestIO(unittest.TestCase):
    """Tests for read/write methods of Structure objects"""

    @classmethod
    def setUpClass(cls):

        cls.pdbpath = TESTDIR / 'data' / 'mini.pdb'
        cls.cifpath = TESTDIR / 'data' / 'mini.cif'
        cls.openmm_structure = app.PDBFile(str(cls.pdbpath))

        # Create temporary folder to write files to
        td = tempfile.mkdtemp()
        cls.tempdir = pathlib.Path(td)

        warnings.simplefilter('ignore')  # blanket warning ignore filter

    @classmethod
    def tearDownClass(cls):
        """Remove tempdir (should be empty)"""

        shutil.rmtree(cls.tempdir.resolve())

    # Reader tests
    def test_default_read(self):
        """Tests reading a structure with default parameters"""

        try:
            s = Structure('dummy', self.openmm_structure)
        except Exception as e:
            self.fail(f'Failed to create Structure from {self.path}: {e}')
        else:
            has_topology = s.topology is not None
            self.assertTrue(has_topology)

            has_positions = s.positions is not None
            self.assertTrue(has_positions)

            has_kdtree = s._kdt is not None
            self.assertTrue(has_kdtree)

    def test_read_noKDTree(self):
        """Checks if flag to NOT build KDTree is respected"""

        s = Structure('dummy', self.openmm_structure, build_kdtree=False)
        no_kdtree = s._kdt is None
        self.assertTrue(no_kdtree)

    # Add
    def test_add_structure(self):
        """Tests adding two structures together"""

        s1 = Structure('A', self.openmm_structure)
        s2 = Structure('B', self.openmm_structure)
        combined = s1 + s2

        size_s1 = len(list(s1.topology.atoms()))
        size_s2 = len(list(s1.topology.atoms()))
        size_combined = len(list(combined.topology.atoms()))
        self.assertEqual(size_s1 + size_s2, size_combined)

    # Copy
    def test_copy_structure(self):
        """Tests copy() method"""

        s = Structure('A', self.openmm_structure)
        sc = s.copy()

        s_at = sorted((a for a in s.topology.atoms()), key=lambda x: x.index)
        sc_at = sorted((a for a in sc.topology.atoms()), key=lambda x: x.index)

        self.assertEqual(s_at, sc_at)
        self.assertIsNone(sc._kdt)

    def test_builtin_copy_structure(self):
        """Tests Python copy.copy() method"""

        from copy import copy

        s = Structure('A', self.openmm_structure)
        sc = copy(s)

        s_at = sorted((a for a in s.topology.atoms()), key=lambda x: x.index)
        sc_at = sorted((a for a in sc.topology.atoms()), key=lambda x: x.index)

        self.assertEqual(s_at, sc_at)
        self.assertIsNone(sc._kdt)

    def test_builtin_deepcopy_structure(self):
        """Tests Python copy.deepcopy() method"""

        from copy import deepcopy

        s = Structure('A', self.openmm_structure)
        sc = deepcopy(s)

        s_at = sorted((a for a in s.topology.atoms()), key=lambda x: x.index)
        sc_at = sorted((a for a in sc.topology.atoms()), key=lambda x: x.index)

        self.assertEqual(s_at, sc_at)
        self.assertIsNone(sc._kdt)

    # Writer tests
    def test_write_PDB(self):
        """Test writing a Structure in PDB format"""

        s = Structure('dummy', self.openmm_structure)
        newfile = self.tempdir / pathlib.Path('test.pdb')
        s.write(str(newfile.resolve()), ftype='pdb')

        # Compare contents
        with newfile.open() as handle:
            new = handle.readlines()

        with open(self.pdbpath) as handle:
            pdb = handle.readlines()

        self.assertEqual(new[1:], pdb[1:])  # remove remark line w/ date.

        # Remove temp file
        newfile.unlink()

    def test_write_mmCIF(self):
        """Test writing a Structure in mmCIF format."""

        s = Structure('dummy', self.openmm_structure)

        newfile = self.tempdir / pathlib.Path('test.cif')
        s.write(str(newfile.resolve()), ftype='cif')

        # Compare contents
        with newfile.open() as handle:
            new = handle.readlines()

        with open(self.cifpath) as handle:
            cif = handle.readlines()

        self.assertEqual(new[2:], cif[2:])  # rm header

        # Remove temp file
        newfile.unlink()

    def test_write_guess_ftype(self):
        """Test writing a Structure without specifying a file type"""

        s = Structure('dummy', self.openmm_structure)

        newfile = self.tempdir / pathlib.Path('test.pdb')
        s.write(str(newfile.resolve()))

        # Read and compare
        with newfile.open() as handle:
            newcontent = handle.readlines()

        with open(self.pdbpath) as original:
            pdbcontent = original.readlines()

        self.assertEqual(newcontent[1:], pdbcontent[1:])

        newfile.unlink()

    def test_write_existing_file(self):
        """Test writing a Structure to an existing file."""

        s = Structure('dummy', self.openmm_structure)

        with self.assertRaises(StructureWriteError):
            newfile = self.tempdir / pathlib.Path('test.pdb')
            newfile.touch()  # create file
            s.write(str(newfile.resolve()))

        newfile.unlink()

    def test_write_error_no_extension(self):
        """Test writing a Structure to file without extension"""

        s = Structure('dummy', self.openmm_structure)

        with self.assertRaises(StructureWriteError):
            newfile = self.tempdir / pathlib.Path('test')
            s.write(str(newfile.resolve()))

    def test_write_error_bad_extension(self):
        """Test writing a Structure to file with unsupported extension"""

        s = Structure('dummy', self.openmm_structure)

        with self.assertRaises(StructureWriteError):
            newfile = self.tempdir / pathlib.Path('test.abc')
            s.write(str(newfile.resolve()))

    def test_write_error_fname_not_string(self):
        """Test writing to a non-string output file name"""

        s = Structure('dummy', self.openmm_structure)

        with self.assertRaises(StructureWriteError):
            dummy_obj = [None, None]
            s.write(dummy_obj)


class TestStructureManipulation(unittest.TestCase):
    """Tests for add/remove methods of Structure objects"""

    @classmethod
    def setUpClass(cls):
        warnings.simplefilter('ignore')  # blanket warning ignore filter

    def test_remove_waters(self):
        """Test removing waters from a PDB file"""

        fpath = TESTDIR / 'data' / 'mini_wHOH.pdb'
        s = ia.read(str(fpath.resolve()))

        n_res_before = len(list(s.topology.residues()))

        try:
            s.remove_waters()
        except Exception as e:
            self.fail(f'Unknown exception when removing waters: {e}')

        n_res_after = len(list(s.topology.residues()))
        n_diff = n_res_before - n_res_after
        self.assertEqual(n_diff, 4)

    def test_remove_unknowns(self):
        """Test removing unknown residues."""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(str(fpath.resolve()))

        n_res_before = s.topology.getNumResidues()

        try:
            s.remove_unknowns()
        except Exception as e:
            self.fail(f'Unknown exception when removing unknowns: {e}')

        n_res_after = s.topology.getNumResidues()
        n_diff = n_res_before - n_res_after
        self.assertEqual(n_diff, 1)  # should remove DRG

    def test_remove_unknowns_exception(self):
        """Test removing unknown residues with exception keyword"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(str(fpath.resolve()))

        n_res_before = s.topology.getNumResidues()

        try:
            s.remove_unknowns(exceptions=['DRG'])
        except Exception as e:
            self.fail(f'Unknown exception when removing unknowns: {e}')

        n_res_after = s.topology.getNumResidues()
        n_diff = n_res_before - n_res_after
        self.assertEqual(n_diff, 0)  # no difference

    def test_add_to_complete_structure(self):
        """Test adding atoms to structure not missing any."""

        fpath = TESTDIR / 'data' / 'mini.pdb'
        s = ia.read(str(fpath.resolve()))

        n_ini_atoms = s.topology.getNumAtoms()

        s.add_missing_heavy_atoms()

        n_end_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_ini_atoms, n_end_atoms)

    def test_addToIncompleteStructure(self):
        """Test adding atoms to structure missing 1 side-chain atom.
        """

        fpath = TESTDIR / 'data' / 'mini_incomplete.pdb'
        s = ia.read(str(fpath.resolve()))

        n_ini_atoms = s.topology.getNumAtoms()

        s.add_missing_heavy_atoms()

        n_end_atoms = s.topology.getNumAtoms()
        self.assertEqual(n_end_atoms - n_ini_atoms, 1)

    # Capping
    def test_add_unknown_n_terminus(self):
        """Test exception when adding unknown capping group (NT)."""

        fpath = TESTDIR / 'data' / 'mini.pdb'
        s = ia.read(str(fpath.resolve()))

        with self.assertRaises(StructureError):
            caps = [('ACE', 'NME'), ('ABC', 'NME')]
            s.add_capping_groups(chain_capping=caps)

    def test_add_unknown_c_terminus(self):
        """Test exception when adding unknown capping group (CT)."""

        fpath = TESTDIR / 'data' / 'mini.pdb'
        s = ia.read(str(fpath.resolve()))

        with self.assertRaises(StructureError):
            caps = [('ACE', 'NME'), ('ACE', 'XYZ')]
            s.add_capping_groups(chain_capping=caps)

    def test_add_default_termini(self):
        """Test adding default termini caps (ACE/NME)."""

        fpath = TESTDIR / 'data' / 'mini.pdb'
        s = ia.read(str(fpath.resolve()))

        n_ini_atoms = s.topology.getNumAtoms()

        s.add_capping_groups()

        # Assert ACE/NME were added
        for chain in s.topology.chains():
            reslist = [r.name for r in chain.residues()]
            n_res, c_res = reslist[0], reslist[-1]
            self.assertEqual(n_res, 'ACE')
            self.assertEqual(c_res, 'NME')

        # Assert method added the right number of atoms
        # 2x ACE + 2x NME = 2x3 + 2x2 = 10
        n_end_atoms = s.topology.getNumAtoms()
        n_diff = n_end_atoms - n_ini_atoms
        self.assertEqual(n_diff, 10)

    def test_add_termini_not_missing_atoms(self):
        """Tests whether add_capping_groups() sticks to termini"""

        fpath = TESTDIR / 'data' / 'mini_incomplete.pdb'
        s = ia.read(str(fpath.resolve()))

        n_ini_atoms = s.topology.getNumAtoms()

        s.add_capping_groups()

        # Explicitly check number of atoms outside ACE/NME
        _caps = set(('ACE', 'NME'))
        noncap_atoms = [1 for r in s.topology.residues()
                        for a in r.atoms() if r.name not in _caps]
        n_noncap_atoms = sum(noncap_atoms)
        self.assertEqual(n_ini_atoms, n_noncap_atoms)

    def test_add_nonmatching_termini(self):
        """Test exception when adding non-matching capping groups."""

        fpath = TESTDIR / 'data' / 'mini.pdb'
        s = ia.read(str(fpath.resolve()))

        with self.assertRaises(StructureError):
            caps = [('ACE', 'NME'), ('ABC', 'NME'), ['ACE', 'NME']]
            s.add_capping_groups(chain_capping=caps)

        with self.assertRaises(StructureError):
            caps = [('ACE', 'NME')]
            s.add_capping_groups(chain_capping=caps)

    def test_add_partial_capping_groups(self):
        """Test adding mixed/partial capping groups"""

        fpath = TESTDIR / 'data' / 'mini.pdb'
        s = ia.read(str(fpath.resolve()))

        n_ini_atoms = s.topology.getNumAtoms()
        chain_res = [[r for r in c.residues()] for c in s.topology.chains()]

        a_nt_res_ini = chain_res[0][0]
        b_ct_res_ini = chain_res[1][-1]

        caps = [(None, 'NME'), ('ACE', None)]
        s.add_capping_groups(chain_capping=caps)

        chain_res = [[r for r in c.residues()] for c in s.topology.chains()]
        a_nt_res = chain_res[0][0]
        a_ct_res = chain_res[0][-1]

        b_nt_res = chain_res[1][0]
        b_ct_res = chain_res[1][-1]

        # Assertions
        self.assertEqual(a_nt_res_ini.name, a_nt_res.name)  # remain same
        self.assertEqual(b_ct_res_ini.name, b_ct_res.name)
        self.assertEqual(a_ct_res.name, 'NME')  # capped
        self.assertEqual(b_nt_res.name, 'ACE')

        c_ter_atoms = {a.name for a in b_ct_res.atoms()}  # charged CT
        self.assertIn('OXT', c_ter_atoms)

        # Number of atoms
        # ACE + NME + C-terminal OXT = 3 + 2 + 1 = 6
        n_end_atoms = s.topology.getNumAtoms()
        n_diff = n_end_atoms - n_ini_atoms
        self.assertEqual(n_diff, 6)

    # Mutate
    def test_simple_mutation(self):
        """Tests mutation method on single residue"""

        fpath = TESTDIR / 'data' / 'mini_noH.pdb'
        s = ia.read(str(fpath.resolve()))

        reslist = [r for r in s.topology.residues()]
        n_ini_atoms = s.topology.getNumAtoms()

        s.mutate(reslist[0], 'ALA')  # ASN-1

        reslist = [r for r in s.topology.residues()]
        n_end_atoms = s.topology.getNumAtoms()
        n_diff = n_ini_atoms - n_end_atoms

        self.assertEqual(reslist[0].name, 'ALA')
        self.assertEqual(n_diff, 3)

    def test_mutate_to_unknown_residue(self):
        """Tests mutation to unknown residue"""

        fpath = TESTDIR / 'data' / 'mini_noH.pdb'
        s = ia.read(str(fpath.resolve()))

        reslist = [r for r in s.topology.residues()]

        with self.assertRaises(StructureError):
            s.mutate(reslist[0], 'FOO')

    def test_mutate_no_adding_atoms(self):
        """Tests whether mutate() does not add other missing atoms"""

        fpath = TESTDIR / 'data' / 'mini_incomplete.pdb'
        s = ia.read(str(fpath.resolve()))

        n_ini_atoms = s.topology.getNumAtoms()
        reslist = [r for r in s.topology.residues()]

        s.mutate(reslist[0], 'ALA')  # ASN-1

        reslist = [r for r in s.topology.residues()]
        n_end_atoms = s.topology.getNumAtoms()
        n_diff = n_ini_atoms - n_end_atoms

        self.assertEqual(reslist[0].name, 'ALA')
        # ASN has 16 atoms, including Hs
        # ALA shares C, N, O, CA, CB (5)
        # mutate does not add Hs
        self.assertEqual(n_diff, 11)

    # KDTree / Search Methods
    def test_get_neighbors_simple(self):
        """Tests simple case for get_neighbors()"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        neighbors = s.get_neighbors(drg)
        self.assertEqual(len(neighbors), 151)

    def test_get_neighbors_atom_list(self):
        """Tests get_neighbors() from list of Atoms"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [a for r in s.topology.residues() for a in r.atoms()
               if r.name == 'DRG']

        neighbors = s.get_neighbors(drg)
        self.assertEqual(len(neighbors), 151)

    def test_get_neighbors_single_atom(self):
        """Tests get_neighbors() from list of Atoms"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [a for r in s.topology.residues() for a in r.atoms()
               if r.name == 'DRG'][0]

        neighbors = s.get_neighbors(drg)
        self.assertEqual(len(neighbors), 25)

    def test_get_neighbors_wrong_entity_type(self):
        """Tests get_neighbors() from weirdo input"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)

        with self.assertRaises(ValueError):
            s.get_neighbors(s)

    def test_get_neighbors_small_radius(self):
        """Tests setting radius for get_neighbors()"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        neighbors = s.get_neighbors(drg, radius=2.5)
        self.assertEqual(len(neighbors), 9)

    def test_get_neighbors_radius_bad(self):
        """Tests get_neighbors() with bad radius value"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        with self.assertRaises(ValueError):
            s.get_neighbors(drg, radius='abc')

    def test_get_neighbors_radius_negative(self):
        """Tests get_neighbors() with negative radius value"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        with self.assertRaises(ValueError):
            s.get_neighbors(drg, radius=-2.5)

    def test_get_neighbors_return_residues(self):
        """Tests get_neighbors() returning residues"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        neighbors = s.get_neighbors(drg, level='residue')
        self.assertEqual(len(neighbors), 21)

    def test_get_neighbors_return_chains(self):
        """Tests get_neighbors() returning chains"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        neighbors = s.get_neighbors(drg, level='chain')
        self.assertEqual(len(neighbors), 1)

    def test_get_neighbors_level_bad(self):
        """Tests get_neighbors() with bad level value"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)
        drg = [r for r in s.topology.residues() if r.name == 'DRG'][0]

        with self.assertRaises(ValueError):
            s.get_neighbors(drg, level='bananas')

    def test_get_neighbors_method_centroid(self):
        """Tests get_neighbors() method as centroid"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)

        drg = [a for r in s.topology.residues() for a in r.atoms()
               if r.name == 'DRG']

        neighbors = s.get_neighbors(drg, method='centroid')
        self.assertEqual(len(neighbors), 16)

    def test_get_neighbors_wrong_method(self):
        """Tests get_neighbors() for wrong method value"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)

        with self.assertRaises(ValueError):
            s.get_neighbors(s, method='bananas')

    def test_get_neighbors_method_reset(self):
        """Tests get_neighbors() method reset for single Atom"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)

        drg = [a for r in s.topology.residues() for a in r.atoms()
               if r.name == 'DRG'][0]

        neighbors = s.get_neighbors(drg, method='centroid')
        self.assertEqual(len(neighbors), 25)  # overrides

    def test_get_neighbor_pairs_simple(self):
        """Tests get_neighbor_pairs() for simple case"""

        fpath = TESTDIR / 'data' / 'toy.pdb'
        s = ia.read(fpath)

        neighbor_pairs = s.get_neighbor_pairs(radius=99.0)
        # 13 atoms = (13*13 - 13) / 2
        self.assertEqual(len(neighbor_pairs), 78)

    def test_get_neighbor_pairs_residue(self):
        """Tests get_neighbor_pairs() returning residue"""

        fpath = TESTDIR / 'data' / 'toy.pdb'
        s = ia.read(fpath)

        neighbor_pairs = s.get_neighbor_pairs(level='residue')
        self.assertEqual(len(neighbor_pairs), 1)  # AB == BA

    def test_get_neighbor_pairs_chain(self):
        """Tests get_neighbor_pairs() returning chain"""

        fpath = TESTDIR / 'data' / 'toy.pdb'
        s = ia.read(fpath)

        neighbor_pairs = s.get_neighbor_pairs(level='chain')
        self.assertEqual(len(neighbor_pairs), 1)  # AB == BA

    def test_get_neighbor_pairs_level_bad(self):
        """Tests get_neighbor_pairs() with bad level value"""

        fpath = TESTDIR / 'data' / 'protein_ligand.pdb'
        s = ia.read(fpath)

        with self.assertRaises(ValueError):
            s.get_neighbor_pairs(level='bananas')


        # class TestReader(unittest.TestCase):
#     """Tests for read() method and associated functions."""

#     def setUp(self):
#         fpath = os.path.join(TESTDIR, 'data', 'protein_ligand.pdb')
#         self.struct = ia.read(fpath)

#     def test_guess_bonds_from_coordinates(self):
#         """Tests __guess_bonds_from_coordinates for DRG ligand."""

#         struct = self.struct
#         # Check DRG residue
#         for r in struct.topology.residues():
#             if r.name == 'DRG':
#                 n_bonds = len(list(drg.bonds()))
#                 self.assertEqual(n_bonds, 26)

#     def test_make_residue_graph(self):
#         """Tests creating a nx.Graph from residue topology.
#         """

#         struct = self.struct
#         reslist = [r for r in struct.topology.residues()]

#         for residue in reslist:
#             g = residue._g

#             n_atoms = len(list(residue.atoms()))
#             n_nodes = g.number_of_nodes()
#             self.assertEqual(n_atoms, n_nodes)

#             n_bonds = len(list(residue.internal_bonds()))
#             n_edges = g.number_of_edges()
#             self.assertEqual(n_bonds, n_edges)

#             # Check C-alpha has 4 edges
#             if residue.name == 'DRG':
#                 continue

#             for idx, atom in enumerate(residue.atoms()):
#                 if atom.name == 'CA':
#                     n_edges = len(g[idx])
#                     self.assertEqual(n_edges, 4)


if __name__ == "__main__":

    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
