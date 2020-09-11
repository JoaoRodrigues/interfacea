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

"""Unit tests for PDB reading module."""

import numpy as np

from interfacea.core.atom import DisorderedAtom
from interfacea.io.pdb import (
    read_pdb,
    PDBParser,
    # PDBParserError,
    # PDBParserWarning
)


# 1ctf: monomer, single-chain, no altlocs or icodes.
def test_reading_1ctf_direct(datadir):
    """Read 1ctf.pdb directly with PDBParser."""

    with (datadir / "1ctf.pdb").open() as handle:
        s = PDBParser(handle)

    # Test Topology
    assert s.topology.natoms == s.topology.nlocs == 554
    assert {a.chain for a in s.topology} == {"A"}
    assert [a.resname for a in s.topology if a.name == "CA"] == [
        "GLU", "PHE", "ASP", "VAL", "ILE", "LEU",
        "LYS", "ALA", "ALA", "GLY", "ALA", "ASN",
        "LYS", "VAL", "ALA", "VAL", "ILE", "LYS",
        "ALA", "VAL", "ARG", "GLY", "ALA", "THR",
        "GLY", "LEU", "GLY", "LEU", "LYS", "GLU",
        "ALA", "LYS", "ASP", "LEU", "VAL", "GLU",
        "SER", "ALA", "PRO", "ALA", "ALA", "LEU",
        "LYS", "GLU", "GLY", "VAL", "SER", "LYS",
        "ASP", "ASP", "ALA", "GLU", "ALA", "LEU",
        "LYS", "LYS", "ALA", "LEU", "GLU", "GLU",
        "ALA", "GLY", "ALA", "GLU", "VAL", "GLU",
        "VAL", "LYS",
    ]

    assert [a.resname for a in s.topology if a.hetatm] == [
        "SO4", "SO4", "SO4", "SO4", "SO4", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH", "HOH", "HOH", "HOH", "HOH", "HOH",
        "HOH",
    ]

    # Test Coordinates
    assert s.xyz.shape == (1, 554, 3)
    assert np.allclose(s.xyz[0, 0, :], [1.274, 13.501, -11.357])
    assert np.allclose(s.xyz.mean(axis=1), [[0.0, 0.0, 0.0]], atol=1e-1)


def test_reading_1ctf_as_path(datadir):
    """Read 1ctf.pdb through read_pdb (as Path)."""

    s = read_pdb(datadir / "1ctf.pdb")

    assert s.topology.natoms == s.topology.nlocs == 554
    assert s.raw_coords.shape == (1, 554, 3)


def test_reading_1ctf_as_str(datadir):
    """Read 1ctf.pdb through read_pdb (as str)."""

    s = read_pdb(str(datadir / "1ctf.pdb"))

    assert s.topology.natoms == s.topology.nlocs == 554
    assert s.raw_coords.shape == (1, 554, 3)


def test_reading_1ctf_as_handle(datadir):
    """Read 1ctf.pdb through read_pdb (as file object)."""

    with (datadir / "1ctf.pdb").open() as handle:
        s = read_pdb(handle)

    assert s.topology.natoms == s.topology.nlocs == 554
    assert s.raw_coords.shape == (1, 554, 3)


def test_reading_1ctf_gzip_as_path(datadir):
    """Read 1ctf.pdb.gz through read_pdb (as path)."""

    s = read_pdb(datadir / "1ctf.pdb.gz")

    assert s.topology.natoms == s.topology.nlocs == 554
    assert s.raw_coords.shape == (1, 554, 3)


def test_reading_1ctf_gzip_as_str(datadir):
    """Read 1ctf.pdb.gz through read_pdb (as str)."""

    s = read_pdb(str(datadir / "1ctf.pdb.gz"))

    assert s.topology.natoms == s.topology.nlocs == 554
    assert s.raw_coords.shape == (1, 554, 3)


def test_reading_1ctf_gzip_as_gzipf(datadir):
    """Read 1ctf.pdb.gz through read_pdb (as gzip file object)."""

    import gzip

    with gzip.open(datadir / "1ctf.pdb.gz", "rt") as handle:
        s = read_pdb(handle)

    assert s.topology.natoms == s.topology.nlocs == 554
    assert s.raw_coords.shape == (1, 554, 3)


# 1brs: multimer, 6 chains, altlocs, no icodes
def test_reading_1brs(datadir):
    """Read 1brs.pdb.gz through read_pdb."""

    s = read_pdb(datadir / "1brs.pdb.gz")

    assert s.topology.natoms == 5151
    assert s.topology.nlocs == 5153
    assert s.raw_coords.shape == (1, 5153, 3)

    chains = {a.chain for a in s.topology}
    assert chains == {"A", "B", "C", "D", "E", "F"}

    atoms = s.topology.atoms
    assert isinstance(atoms[2686], DisorderedAtom)
    assert len(atoms[2686]) == 2  # number of altlocs


# 1brs: NMR structure, multimer, 2 chains
def test_reading_1ggr(datadir):
    """Read 1ggr.pdb.gz through read_pdb."""

    s = read_pdb(datadir / "1ggr.pdb.gz")

    assert s.topology.natoms == 3595
    assert s.topology.nlocs == 3595
    assert s.raw_coords.shape == (3, 3595, 3)

    chains = {a.chain for a in s.topology}
    assert chains == {"A", "B"}

    assert np.allclose(s.raw_coords[0, 0, :], [51.382, 49.015, 12.266])
    assert np.allclose(s.raw_coords[1, 0, :], [51.340, 48.946, 12.186])  # yes, same.
    assert np.allclose(s.raw_coords[2, 0, :], [51.340, 48.946, 12.186])  # yes, same


def test_reading_1ggr_first_model(datadir):
    """Read only first model of 1ggr.pdb.gz through read_pdb."""

    s = read_pdb(datadir / "1ggr.pdb.gz", nmodels=1)

    assert s.topology.natoms == 3595
    assert s.topology.nlocs == 3595
    assert s.raw_coords.shape == (1, 3595, 3)

    chains = {a.chain for a in s.topology}
    assert chains == {"A", "B"}

    assert np.allclose(s.raw_coords[0, 0, :], [51.382, 49.015, 12.266])
