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

import gzip
import io
import pathlib

import pytest

from interfacea.core.atom import Atom
from interfacea.io.pdb import (
    PDBParser,
    PDBParserError,
    PDBParserWarning
)


datadir = pathlib.Path('tests') / 'data'


def test_no_arguments():
    """Throws error when creating PDBParser without arguments"""

    with pytest.raises(PDBParserError):
        _ = PDBParser()


def test_wrong_arguments():
    """Throws error when creating PDBParser with wrong arguments"""

    with pytest.raises(PDBParserError):
        _ = PDBParser(file='something', pdbid='else')


@pytest.mark.parametrize(
    'datafile',
    (
        datadir / 'tripeptide.pdb',
        datadir / '1brs.pdb',
        datadir / '1ggr.pdb',
    )
)
def test_load_local_files(datafile):
    """Create a new PDBParser from a file on disk"""

    p = PDBParser(datafile)
    assert isinstance(p.source, io.TextIOWrapper)


@pytest.mark.usefixtures('has_internet_connectivity')
def test_load_remote_file():
    """Create a new PDBParser from a file on disk"""

    p = PDBParser(pdbid='1brs')
    assert isinstance(p.source, gzip.GzipFile)


@pytest.mark.parametrize(
    'datafile',
    (
        datadir / 'tripeptide.pdb',
        datadir / '1ggr.pdb',
        datadir / '1brs.pdb'
    )
)
def test_assert_types(datafile):
    """Asserts the types of the attributes of the abstract class."""

    p = PDBParser(datafile)
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


@pytest.mark.parametrize(
    ('datafile', 'natoms', 'nmodels'),
    (
        (datadir / 'tripeptide.pdb', 107, 1),
        (datadir / '1ggr.pdb', 3595, 3),
        (datadir / '1brs.pdb', 5153, 1)
    )
)
def test_parse(datafile, natoms, nmodels):
    """Successfully parses test PDB files."""

    p = PDBParser(datafile)
    p.parse()

    assert len(p.atoms) == len(p.coords) == natoms
    assert sum(map(len, p.coords)) == natoms * nmodels
    assert len(p.coords[0]) == nmodels


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

    p = PDBParser(datafile)
    with pytest.raises(PDBParserError, match=emsg):
        p.parse()


def test_auto_assign_elements():
    """Assigns elements correctly."""

    p = PDBParser(datadir / 'missingelement.pdb')

    with pytest.warns(PDBParserWarning, match='Could not assign element'):
        p.parse()

    elems = [a.element.symbol for a in p.atoms]
    assert elems == ['N', 'H', 'H', 'H', 'C', 'X', 'P', 'X']


@pytest.mark.parametrize(
    ('datafile', 'n_hetatm'),
    (
        (datadir / 'tripeptide.pdb', 0),
        (datadir / '1ggr.pdb', 4),
        (datadir / '1brs.pdb', 513)
    )
)
def test_hetatm(datafile, n_hetatm):
    """Identifies and assigns HETATM correctly."""

    p = PDBParser(datafile)
    p.parse()

    hetatm = list(filter(lambda a: a.hetatm, p.atoms))
    assert len(hetatm) == n_hetatm
