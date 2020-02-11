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
Unit tests for core.structure classes:
    - Atom
    - DisorderedAtom
    - Structure
"""

import pathlib

import numpy as np
import pytest

from interfacea.io import read as read_structure
from interfacea.io.atomrecord import AtomRecord
from interfacea.core.structure import (
    Atom,
    DisorderedAtom,
    Structure
)
import interfacea.exceptions as e


def test_dummy_structure():
    """Successfully creates a small structure from Atoms"""
    _ = Structure(
        'dummy',
        [Atom('N', 0)],
        np.array([[[1.0, 1.0, 1.0]]]),
    )


# Structure Class
class TestStructure:

    def setup(self):

        rootdir = pathlib.Path(".")
        datadir = rootdir / "tests" / "data"

        self.defaultmdl = str(datadir / 'pdb' / 'default.pdb')
        self.multimodel = str(datadir / 'pdb' / 'multimodel.pdb')

    # class Atom
    @pytest.mark.parametrize(
        "attr_dict",
        (
            # Mandatory Args.
            {
                'name': 'CA',
                'serial': 0
            },
            # Optional Parameters
            {
                'name': 'CA',
                'serial': 0,
                'resid': 140,
                'chain': 'A',
            },
        )
    )
    def test_Atom_optional_args(self, attr_dict):
        """Manually creates an instance of an Atom object"""
        a = Atom(**attr_dict)
        for key, value in attr_dict.items():
            assert getattr(a, key) == value

    @pytest.mark.parametrize(
        "attrs",
        (
            {
                'serial': 0,
                'model': 1,
                'rectype': 'ATOM',
                'name': 'CA',
                'altloc': '',
                'resname': 'ALA',
                'chain': 'A',
                'resid': 15,
                'icode': '',
                'x': 1.0,
                'y': 1.0,
                'z': 1.0,
                'occ': 1.0,
                'b': 0.0,
            },
        )
    )
    def test_Atom_from_AtomRecord(self, attrs):
        """Manually creates an instance of an Atom object"""
        atomrecord = AtomRecord(*attrs.values())
        a = Atom.from_atomrecord(atomrecord)
        for key, value in attrs.items():
            if key in ('x', 'y', 'z'):
                assert getattr(a, key, None) is None
            else:
                assert getattr(a, key) == value

    @pytest.mark.parametrize(
        "attr_dict",
        (
            {'name': 'CA'},
            {'serial': 0},
        )
    )
    def test_Atom_wrong_args(self, attr_dict):
        """Raises error when Atom() misses mandatory arguments"""
        with pytest.raises(TypeError):
            _ = Atom(**attr_dict)

    def test_Atom_str(self):
        """Successfully outputs Atom string representation"""
        a = Atom('N', 0)
        assert a.__str__() == "<Atom name=N serial=0>"

    def test_Atom_unbound_parent(self):
        """Unbound Atom has parent None"""
        a = Atom('N', 0)
        assert a.parent is None

    def test_Atom_unbound_coords(self):
        """Unbound Atom raises error on coord access"""
        a = Atom('N', 0)
        with pytest.raises(AttributeError):
            _ = a.coords

    # class DisorderedAtom
    def test_DisorderedAtom(self):
        """Successfully creates DisorderedAtom"""
        aa = Atom('N', 0, altloc='A', occ=0.8)
        ab = Atom('C', 1, altloc='B', occ=0.2)

        atoms = [aa, ab]

        da = DisorderedAtom()

        with pytest.raises(AttributeError):
            _ = da.name  # no children yet

        da.from_list(atoms)

        assert da.nlocs == 2
        assert da.selected_child == aa
        assert da.name == 'N' and da.serial == 0

        with pytest.raises(AttributeError):
            _ = da.missing_attr  # unknown

        # Change attributes (forwarding of setattr)
        da.name = 'O'
        assert da.name == 'O'

    def test_DisorderedAtom_add(self):
        """Successfully adds single Atoms to DisorderedAtom """
        aa = Atom('N', 0, altloc='A', occ=0.2)
        ab = Atom('C', 1, altloc='B', occ=0.8)

        da = DisorderedAtom()

        da.add(aa)
        assert da.nlocs == 1
        assert da.name == 'N'

        da.add(ab)
        assert da.nlocs == 2
        assert da.name == 'C'  # occ(b) > occ(a)

        with pytest.raises(e.DuplicateAltLocError):
            da.add(aa)

    def test_DisorderedAtom_select(self):
        """Successfully swaps DisorderedAtom locs"""
        aa = Atom('N', 0, altloc='A', occ=0.8)
        ab = Atom('C', 1, altloc='B', occ=0.2)

        atoms = [aa, ab]

        da = DisorderedAtom()
        da.from_list(atoms)

        assert da.selected_child == aa
        assert da.name == 'N' and da.serial == 0

        da.select('B')

        assert da.selected_child == ab
        assert da.name == 'C' and da.serial == 1

        with pytest.raises(KeyError):
            da.select('Z')

    # class Structure
    def test_Structure_io_read(self):
        """io.read() successfully builds Structure"""

        s = read_structure(self.defaultmdl)  # io.read()

        assert s.name == 'default.pdb'
        assert s.natoms == 107
        assert s.nmodels == 1
        assert s.precision == np.float16

        # Check DisorderedAtom
        # Check atoms are bound
        for atom in s.atoms:
            if atom.serial == 18:
                assert isinstance(atom, DisorderedAtom)
                assert atom.nlocs == 2
            else:
                assert isinstance(atom, Atom)

            assert atom.parent is s

    def test_Structure_io_read_noaltlocs(self):
        """io.read() successfully builds Structure (no altloc)"""

        s = read_structure(self.defaultmdl, discard_altloc=True)

        assert s.name == 'default.pdb'
        assert s.natoms == 106
        assert s.nmodels == 1
        assert s.precision == np.float16

        # Check no DisorderedAtom and best occ
        for atom in s.atoms:
            if atom.serial == 18:
                assert atom.occ == 0.8
            assert isinstance(atom, Atom)
            assert atom.parent is s

    def test_Structure_multimodel(self):
        """io.read() successfully builds multi-model Structure"""

        s = read_structure(self.multimodel)

        assert s.name == 'multimodel.pdb'
        assert s.natoms == 106
        assert s.nmodels == 2
        assert s.precision == np.float16

    def test_Structure_multimodel_selection(self):
        """Structure.select successfully changes active model"""

        s = read_structure(self.multimodel)

        # Test model swapping
        assert s.model == 0
        assert s.atoms[0].coords[0] < 100.0

        s.model = 1

        assert s.model == 1
        assert s.atoms[0].coords[0] > 100.0

        with pytest.raises(ValueError):
            s.model = -1

    def test_Structure_access_coords(self):
        """Structure.coords methods return adequate arrays"""

        s = read_structure(self.multimodel)

        assert s.coords.shape == (106, 3)
        assert s.full_coords.shape == (2, 106, 3)

        # Setter raises error
        with pytest.raises(NotImplementedError):
            s.coords = [1.0, 1.0, 1.0]

        with pytest.raises(NotImplementedError):
            s.full_coords = [1.0, 1.0, 1.0]
