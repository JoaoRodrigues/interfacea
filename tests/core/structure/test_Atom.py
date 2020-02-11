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
"""

import pytest

from interfacea.core.structure import (
    Atom,
)
from interfacea.io.atomrecord import AtomRecord


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
def test_Atom_optional_args(attr_dict):
    """Creates instance of Atom manually."""
    a = Atom(**attr_dict)
    for key, value in attr_dict.items():
        assert getattr(a, key) == value


def test_Atom_from_AtomRecord():
    """Creates instance of Atom from AtomRecord."""

    attrs = {
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
    }

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
def test_Atom_wrong_args(attr_dict):
    """Raises error when instantiating Atom without mandatory arguments"""
    with pytest.raises(TypeError):
        _ = Atom(**attr_dict)


def test_Atom_unbound_parent():
    """Returns None for Atom.parent when unbound."""
    a = Atom('N', 0)
    assert a.parent is None


def test_Atom_unbound_coords():
    """Raises error when accessing Atom.coords when unbound."""
    a = Atom('N', 0)
    with pytest.raises(AttributeError):
        _ = a.coords
