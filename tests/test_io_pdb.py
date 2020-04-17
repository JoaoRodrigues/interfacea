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
Unit tests for PDBParser.
"""

import collections
import pathlib

import pytest

# For type checking
from interfacea.chemistry.elements import Element
from interfacea.core.atom import Atom

from interfacea.io.pdb import (
    PDBParser,
    PDBParserError,
    PDBParserWarning
)


datadir = pathlib.Path('tests') / 'data'


# Setup Fixture
@pytest.fixture(
    scope='module',
    params=[
        (datadir / 'tripeptide.pdb', 107, 1, 0),
        (datadir / '1ggr.pdb', 3595, 3, 4),
        (datadir / '1brs.pdb', 5153, 1, 513)
    ])
def pdbfile(request):
    s = collections.namedtuple(
        'Structure',
        ['path', 'natoms', 'nmodels', 'n_hetatm']
    )
    return s(*request.param)


def test_no_arguments():
    """Throws error when creating PDBParser without arguments"""

    with pytest.raises(TypeError):
        _ = PDBParser()


def test_wrong_argtype():
    """Throws error when creating PDBParser with wrong argument type"""

    with pytest.raises(TypeError):
        _ = PDBParser(file=12345)


def test_load_text_files(pdbfile):
    """Create a new PDBParser from local files"""

    with pdbfile.path.open('rt') as handle:
        _ = PDBParser(handle)


def test_assert_types(pdbfile):
    """Asserts the types of the attributes of the abstract class."""

    with pdbfile.path.open('rt') as handle:
        p = PDBParser(handle)
        p.parse()

    atoms = p.atoms
    assert isinstance(atoms, list)
    assert isinstance(atoms[0], Atom)
    assert len(set(map(type, atoms))) == 1

    coords = p.coords
    assert isinstance(coords, list)  # 1 per atom
    assert isinstance(coords[0], list)  # 1 per model
    assert isinstance(coords[0][0], list)  # x,y,z
    assert len(coords[0][0]) == 3
    assert isinstance(coords[0][0][0], float)  # x
    assert len(set(map(type, coords[0]))) == 1


def test_parse(pdbfile):
    """Successfully parses test PDB files."""

    with pdbfile.path.open('rt') as handle:
        p = PDBParser(handle)
        p.parse()

    assert len(p.atoms) == len(p.coords) == pdbfile.natoms
    assert sum(map(len, p.coords)) == pdbfile.natoms * pdbfile.nmodels
    assert len(p.coords[0]) == pdbfile.nmodels


@pytest.mark.parametrize(
    ('datafile', 'emsg'),
    (
        (datadir / 'badatom.pdb', 'line 3'),
        (datadir / 'badmodel.pdb', 'Unexpected MODEL record'),
        (datadir / 'badendmdl.pdb', 'Unexpected ENDMDL record'),
        (datadir / 'duplicatedatom.pdb', 'is duplicated on line'),

    )
)
def test_parse_errors(datafile, emsg):
    """Successfully parses test PDB files."""

    with pytest.raises(PDBParserError, match=emsg):
        with datafile.open('rt') as handle:
            p = PDBParser(handle)
            p.parse()


def test_auto_assign_elements():
    """Assigns elements correctly."""

    fpath = datadir / 'missingelement.pdb'
    with pytest.warns(PDBParserWarning, match='Could not assign element'):
        with fpath.open('rt') as handle:
            p = PDBParser(handle)
            p.parse()

    elem_list = ['N', 'H', 'H', 'H', 'C', 'X', 'P', 'X']
    for i, a in enumerate(p.atoms):
        assert isinstance(a.element, Element)
        assert a.element.symbol == elem_list[i]


def test_hetatm(pdbfile):
    """Identifies and assigns HETATM correctly."""

    with pdbfile.path.open('rt') as handle:
        p = PDBParser(handle)
        p.parse()

    hetatm = list(filter(lambda a: a.hetatm, p.atoms))
    assert len(hetatm) == pdbfile.n_hetatm
