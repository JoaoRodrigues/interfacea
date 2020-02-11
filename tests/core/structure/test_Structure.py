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
    - Structure
"""

import numpy as np
import pytest

from interfacea.io import read as io_read
from interfacea.core.structure import (
    Atom,
    DisorderedAtom,
    Structure
)

from . import datadir


defaultpdb = str(datadir / 'pdb' / 'default.pdb')
multimodel = str(datadir / 'pdb' / 'multimodel.pdb')


# Dummy Structure
def test_create_Structure():
    """Creates instance of Structure from Atoms and np.array"""

    aa = Atom('N', 0)

    s = Structure(
        'dummy',
        [aa],
        np.array([[[1.0, 1.0, 1.0]]]),
    )

    assert s.natoms == 1
    assert s.atoms[0] == aa


# class Structure
def test_create_Structure_from_file():
    """Creates instance of Structure from file on disk (via io.read)"""

    s = io_read(defaultpdb)

    assert s.name == 'default.pdb'
    assert s.natoms == 107
    assert s.nmodels == 1
    assert s.precision == np.float16

    assert s.coords.shape == (107, 3)
    assert s.full_coords.shape == (1, 107, 3)


def test_create_Structure_withDisorder():
    """Creates DisorderedAtoms when building Structure"""

    s = io_read(defaultpdb)

    assert s.natoms == 107

    for atom in s.atoms:
        if atom.serial == 18:
            assert isinstance(atom, DisorderedAtom)
            assert atom.nlocs == 2
        else:
            assert isinstance(atom, Atom)

        assert atom.parent is s  # is bound?


def test_create_Structure_noDisorder():
    """Ignore multiple altlocs when building Structure"""

    s = io_read(defaultpdb, discard_altloc=True)

    assert s.natoms == 106

    for atom in s.atoms:
        if atom.serial == 18:
            assert atom.occ == 0.8
        assert isinstance(atom, Atom)


def test_create_Structure_ensemble():
    """Create instance of Structure with multiple models"""

    s = io_read(multimodel)

    assert s.name == 'multimodel.pdb'
    assert s.natoms == 106
    assert s.nmodels == 2
    assert s.precision == np.float16

    assert s.coords.shape == (106, 3)
    assert s.full_coords.shape == (2, 106, 3)


def test_Structure_get_atom():
    """Returns Atom object on Structure.get_atom"""

    s = io_read(defaultpdb)

    with pytest.raises(IndexError):
        s.atoms[106]  # does not exist

    atom = s.get_atom(106)
    assert atom.serial == 106

    # Return same object: DisorderedAtom
    assert s.get_atom(18) == s.get_atom(19)


def test_Structure_select_model():
    """Changes coordinate set on Structure.model selection"""

    s = io_read(multimodel)

    assert s.model == 0
    assert s.atoms[0].coords[0] < 100.0

    s.model = 1  # swap

    assert s.model == 1
    assert s.atoms[0].coords[0] > 100.0


@pytest.mark.parametrize(
    "model_no",
    (-1, 10)
)
def test_Structure_select_model_fail(model_no):
    """Raises error on invalid model number when selecting"""

    s = io_read(multimodel)
    with pytest.raises(ValueError):
        s.model = model_no


def test_Structure_coords_setter_fail():
    """Raises error when setting (full) coordinates"""

    s = io_read(multimodel)
    with pytest.raises(NotImplementedError):
        s.coords = [1.0, 1.0, 1.0]

    with pytest.raises(NotImplementedError):
        s.full_coords = [1.0, 1.0, 1.0]
