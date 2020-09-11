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

"""Unit tests for core Structure object."""

import copy
import gzip

import numpy as np
import pytest

from interfacea.core.structure import Structure, StructureError
from interfacea.io.pdb import (
    PDBParser,
)


@pytest.fixture(scope="module")
def pdb1ggr(datadir):
    """Fixture to load structure."""

    fpath = datadir / "1ggr.pdb.gz"

    with gzip.open(fpath, "rt") as handle:
        return PDBParser(handle)


def test_load_structure(pdb1ggr):
    """Create structure from parsed file."""

    parser = pdb1ggr
    s = Structure(topology=parser.topology, xyzdata=parser.xyz)
    assert s.raw_coords.shape == (3, 3595, 3)


def test_select_model(pdb1ggr):
    """Select between models in a structure."""

    parser = pdb1ggr
    s = Structure(topology=parser.topology, xyzdata=parser.xyz)

    assert s.model == 0
    assert np.allclose(s.coords[0, :], [51.382, 49.015, 12.266])

    s.model = 1
    assert s.model == 1
    assert np.allclose(s.coords[0, :], [51.340, 48.946, 12.186])

    s.model = 2
    assert s.model == 2
    assert np.allclose(s.coords[0, :], [51.340, 48.946, 12.186])  # not a typo

    with pytest.raises(StructureError) as excinfo:
        s.model = 3
    assert "Model index must be between 0 and 2" in str(excinfo)


def test_bad_topology(pdb1ggr):
    """Raise error when Topology has inconsistent indexes."""

    parser = pdb1ggr
    top = copy.deepcopy(parser.topology)
    for a in top:
        a.index = 1

    with pytest.raises(StructureError):
        Structure(topology=top, xyzdata=parser.xyz)


def test_bad_coordinates(pdb1ggr):
    """Raise error when coordinates are incorrectly formatted."""

    parser = pdb1ggr

    _ = Structure(topology=parser.topology, xyzdata=parser.xyz.tolist())

    with pytest.raises(StructureError) as excinfo:
        xyz = np.array([])
        Structure(topology=parser.topology, xyzdata=xyz)
    assert "Coordinate array has wrong dimensions" in str(excinfo)

    with pytest.raises(StructureError) as excinfo:
        Structure(topology=parser.topology, xyzdata=parser.xyz[:, :, :2])
    assert "Wrong coordinate array shape" in str(excinfo)

    with pytest.raises(StructureError) as excinfo:
        xyz = np.array([])
        Structure(topology=parser.topology, xyzdata=parser.xyz[:, :10, :])
    assert "Topology size does not match coordinates" in str(excinfo)
