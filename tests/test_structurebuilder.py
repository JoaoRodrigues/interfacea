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
Unit tests for StructureBuilder.
"""

import pytest

from interfacea.core.atom import DisorderedAtom, Atom
from interfacea.io.structurebuilder import StructureBuilder


@pytest.fixture(scope='module')
def atomdata():
    """Sets up a list entries mimicking parsed atom data."""

    # tuple(name, coords, atomdata)
    return [
        ('C', [0.0, 0.0, 0.0], {'resid': 1, 'chain': 'A', 'icode': ''}),
        ('C', [1.0, 0.0, 0.0], {'resid': 1, 'chain': 'A', 'icode': ''}),
        ('N', [0.0, 1.0, 0.0], {'resid': 1, 'chain': 'A', 'icode': ''}),
        ('O', [0.0, 0.0, 1.0], {'resid': 1, 'chain': 'A', 'icode': ''})
    ]


def test_create_error():
    """Throws error when creating StructureBuilder without name"""

    with pytest.raises(TypeError):
        _ = StructureBuilder()


def test_create():
    """Create a new StructureBuilder"""

    sb = StructureBuilder(name='Test')
    assert sb.name == 'Test'
    assert not sb.atoms
    assert not sb.coord
    assert not sb.params  # exists but is empty


def test_create_wparams():
    """Create a new StructureBuilder with params"""

    sb = StructureBuilder(name='Test', skip_altloc=False)
    assert sb.name == 'Test'
    assert not sb.atoms
    assert not sb.coord
    assert sb.params == {'skip_altloc': False}


def test_add_single_atom(atomdata):
    """Add an Atom to StructureBuilder"""

    sb = StructureBuilder(name='mystructure')
    atom = Atom(name='C')
    sb._add_single_atom(atom)
    assert len(sb.atoms) == 1


def test_add_disordered_atom():
    """Makes DisorderedAtom from Atom in StructureBuilder"""

    sb = StructureBuilder(name='mystructure')
    a = Atom(name='C')
    sb.atoms = [a]
    sb._add_disordered_atom(0, Atom(name='C'))
    assert len(sb.atoms) == 1


def test_add_to_disordered_atom():
    """Add altlocs to DisorderedAtom in StructureBuilder"""

    sb = StructureBuilder(name='mystructure')
    da = DisorderedAtom()
    da.add(Atom(name='C'))
    sb.atoms = [da]

    sb._add_disordered_atom(0, Atom(name='C'))
    assert len(sb.atoms) == 1


def test_add_atom(atomdata):
    """Adds a series of atoms to StructureBuilder (default)"""

    sb = StructureBuilder(name='mystructure')
    for name, _, metadata in atomdata:
        sb.add_atom(name, metadata)

    assert len(sb.atoms) == 3
    assert(sum(1 for a in sb.atoms if isinstance(a, DisorderedAtom)) == 1)


def test_add_atom_skipaltloc(atomdata):
    """Adds a series of atoms to StructureBuilder (skip_altloc)"""

    sb = StructureBuilder(name='mystructure', skip_altloc=True)
    for name, _, metadata in atomdata:
        sb.add_atom(name, metadata)

    assert len(sb.atoms) == 3
    assert(sum(1 for a in sb.atoms if isinstance(a, DisorderedAtom)) == 0)


def test_clean(atomdata):
    """Adds a series of atoms to StructureBuilder (skip_altloc)"""

    sb = StructureBuilder(name='mystructure')
    for name, _, metadata in atomdata:
        sb.add_atom(name, metadata)

    sb.clear()
    assert not sb.atoms
    assert not sb.coord
    assert not sb.params
