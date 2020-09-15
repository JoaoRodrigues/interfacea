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

"""Unit tests for KDTree wrapper code and methods."""

import numpy as np
import pytest

from interfacea.io import read
from interfacea.geometry.kdtree import KDTree


@pytest.fixture(scope="module")
def pdb1ggr(datadir):
    """Fixture to load structure."""

    fpath = datadir / "1ggr.pdb.gz"
    return read(fpath)


def test_create_kdtree_instance(pdb1ggr):
    """Create a KDTree instance from a Structure."""

    _ = KDTree(pdb1ggr)


def test_kdtree_find_neighbors_atom(pdb1ggr):
    """Test find_neighbors() with Atom."""

    kdt = KDTree(pdb1ggr)

    atom = pdb1ggr.topology.atoms[48]
    nb = kdt.find_neighbors(atom, radius=2.0)

    # pymol: 1ggr_0001 w. 2A of (n. CA and i. 22 and c. A)
    # 5 neighbors including self.
    assert len(nb) == 4
    assert {(n.name, n.resid, n.chain) for n in nb} == \
        {("N", 22, "A"), ("C", 22, "A"), ("CB", 22, "A"), ("HA", 22, "A")}


def test_kdtree_find_neighbors_array(pdb1ggr):
    """Test find_neighbors() with np.ndarray."""

    pt = np.array([47.531, 44.605, 21.626])  # CA of ILE22 of chain A

    kdt = KDTree(pdb1ggr)
    nb = kdt.find_neighbors(pt, radius=2.0)

    assert len(nb) == 5
    assert {(n.name, n.resid) for n in nb} == \
        {("CA", 22), ("N", 22), ("C", 22), ("CB", 22), ("HA", 22)}


def test_kdtree_find_neighbors_sort(pdb1ggr):
    """Test find_neighbors() with sort=True."""

    pt = np.array([47.531, 44.605, 21.626])  # CA of ILE22 of chain A

    kdt = KDTree(pdb1ggr)
    nb = kdt.find_neighbors(pt, radius=2.0, sort=True)

    assert len(nb) == 5
    assert [(n.name, n.resid) for n in nb] == \
        [("CA", 22), ("HA", 22), ("N", 22), ("C", 22), ("CB", 22)]


def test_kdtree_find_neighbors_distances(pdb1ggr):
    """Test find_neighbors() with distances=True."""

    pt = np.array([47.531, 44.605, 21.626])  # CA of ILE22 of chain A

    kdt = KDTree(pdb1ggr)
    nb = kdt.find_neighbors(pt, radius=2.0, distances=True)

    expected = {
        "CA": 0.00,
        "HA": 1.08,
        "N": 1.46,
        "C": 1.51,
        "CB": 1.55
    }

    assert len(nb) == 5
    for (a, d) in nb:
        assert round(d, 2) == expected[a.name]


def test_fail_bucket_size(pdb1ggr):
    """Raise error on invalid bucket_size parameter."""

    with pytest.raises(ValueError):
        KDTree(pdb1ggr, bucket_size=0)


def test_kdtree_fail_radius(pdb1ggr):
    """Raise error on invalid radius parameter."""

    kdt = KDTree(pdb1ggr)

    atom = pdb1ggr.topology.atoms[48]
    with pytest.raises(ValueError):
        kdt.find_neighbors(atom, radius=-1.0)


def test_kdtree_fail_center(pdb1ggr):
    """Raise error on invalid center argument."""

    kdt = KDTree(pdb1ggr)

    with pytest.raises(TypeError):
        kdt.find_neighbors("1.0, 2.0, 3.0", radius=1.0)

    pt = np.array([47.531, 44.605, 21.626, 1.234])
    with pytest.raises(ValueError):
        kdt.find_neighbors(pt, radius=1.0)
