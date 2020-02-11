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
Unit tests for KDTree implementation via Structure/Atom.
"""

import pytest

from interfacea.io import read as io_read

from . import datadir

defaultpdb = str(datadir / 'pdb' / 'default.pdb')


# KDTree Tests
def test_Structure_create_kdtree():
    """Creates a KDTree when building Structure"""

    s = io_read(defaultpdb)

    assert s._kdtree is not None
    assert hasattr(s._kdtree, 'search')


def test_Atom_neighbor_search():
    """Returns neighbors of Atom objects"""

    s = io_read(defaultpdb)

    center = s.get_atom(0)
    nlist = list(center.neighbors(radius=2.0))

    assert len(nlist) == 4
    assert sorted([n.name for n, _ in nlist]) == \
        ['CA', 'H', 'H2', 'H3']


@pytest.mark.parametrize(
    "r",
    (-1, 0)
)
def test_Atom_neighbor_search_fail(r):
    """Fails when radius in Atom.neighbors is zero or negative"""

    s = io_read(defaultpdb)

    center = s.get_atom(0)
    with pytest.raises(ValueError):
        _ = list(center.neighbors(radius=r))  # generator
