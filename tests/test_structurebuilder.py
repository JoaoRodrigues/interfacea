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

from interfacea.core.atom import DisorderedAtom
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

    sb = StructureBuilder(name='Test', use_altloc=False)
    assert sb.name == 'Test'
    assert not sb.atoms
    assert not sb.coord
    assert sb.params == {'use_altloc': False}


def test_add_atom(atomdata):
    """Add Atoms to new StructureBuilder"""

    sb = StructureBuilder(name='mystructure')
    for name, coords, metadata in atomdata:
        sb.add_atom(name, coords, metadata)

    n_disordered = sum(isinstance(a, DisorderedAtom) for a in sb.atoms)
    assert n_disordered == 1
    assert len(sb.atoms) == 3  # != from serial for DisorderedAtoms


def test_add_atom_skipaltloc(atomdata):
    """Add Atoms to new StructureBuilder (skip_altloc=True)"""

    sb = StructureBuilder(name='mystructure', skip_altloc=True)
    for name, coords, metadata in atomdata:
        sb.add_atom(name, coords, metadata)

    n_disordered = sum(isinstance(a, DisorderedAtom) for a in sb.atoms)
    assert n_disordered == 0
    assert len(sb.atoms) == 3
