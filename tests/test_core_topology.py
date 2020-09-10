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

"""Unit tests for Topology class."""

import copy

import networkx as nx
import pytest

from interfacea.core.atom import Atom, DisorderedAtom
from interfacea.core.topology import Topology, TopologyError


def test_instantiate_success():
    """Instantiate the Topology class."""

    t = Topology()

    assert t.atoms == []
    assert t.natoms == 0
    assert t.nlocs == 0


def test_str_dunder():
    """Test str(Topology)."""

    t = Topology()
    assert t.__str__() == "Topology with 0 atoms"


# Copy
def test_copy_shallow():
    """Make a shallow copy of the topology."""

    atoms = [Atom(name="CA"), Atom(name="CB"), Atom(name="CG"), Atom(name="CD")]
    bonds = nx.Graph([(0, 1), (1, 2), (2, 3)])  # 4 nodes

    t = Topology()
    for a in atoms:
        t.add_atom(a)
    t.bonds = bonds

    t2 = copy.copy(t)
    assert t is not t2
    assert t.bonds is t2.bonds
    for a1, a2 in zip(t.atoms, t2.atoms):
        assert a1 is a2


def test_copy_deep():
    """Deep copy of the topology."""

    atoms = [Atom(name="CA"), Atom(name="CB"), Atom(name="CG"), Atom(name="CD")]
    bonds = nx.Graph([(0, 1), (1, 2), (2, 3)])  # 4 nodes

    t = Topology()
    for a in atoms:
        t.add_atom(a)
    t.bonds = bonds

    t2 = copy.deepcopy(t)
    assert t is not t2

    assert t.bonds is not t2.bonds
    assert nx.is_isomorphic(t.bonds, t2.bonds) is True

    for a1, a2 in zip(t.atoms, t2.atoms):
        assert a1 is not a2
        assert a1 == a2


# Adding atoms
def test_add_atoms():
    """Add Atoms to Topology."""

    a = Atom(name="CA")

    t = Topology()
    t.add_atom(a)

    assert t.natoms == 1
    assert t.atoms[0] == a
    assert a.index == 0

    a = Atom(name="CB")
    t.add_atom(a)

    assert t.natoms == 2
    assert t.nlocs == 2
    assert t.atoms[1] == a
    assert a.index == 1


def test_add_disorderedatom():
    """Add DisorderedAtom to Topology."""

    da = DisorderedAtom()
    a = Atom(name="CA")
    da.add(a)

    t = Topology()
    t.add_atom(da)

    assert t.natoms == 1
    assert t.nlocs == 1
    assert t.atoms[0] == da
    assert a.index == 0


def test_make_disorderedatom_from_atom():
    """Make DisorderedAtom when adding duplicate Atoms to Topology."""

    a = Atom(name="CA", altloc="A")
    b = Atom(name="CA", altloc="B")

    t = Topology()
    t.add_atom(a)
    t.add_atom(b)

    assert t.natoms == 1  # bundled in one DisorderedAtom
    assert t.nlocs == 2


def test_add_atom_ignore_altloc(caplog):
    """Ignore altlocs when adding duplicate Atoms to Topology."""

    a = Atom(name="CA", altloc="A")

    t = Topology()
    t.add_atom(a)
    t.add_atom(a, disorder=False)

    assert t.natoms == 1
    assert t.nlocs == 1

    assert len(caplog.records) == 1
    for record in caplog.records:
        assert record.levelname == "WARNING"
        assert "Ignored alternate location of atom" in caplog.text


# Get/Set methods
def test_getset_atoms():
    """Retrieve/Set an Atom from/in the Topology."""

    a1 = Atom(name="CA")
    a2 = Atom(name="CB")

    t = Topology()
    t.add_atom(a1)

    assert t[0] is a1
    t[0] = a2
    assert t[0] is a2

    with pytest.raises(IndexError):
        t[2]

    with pytest.raises(IndexError):
        t[2] = a2


# Iteration
def test_iter():
    """Iterate over Topology atoms."""

    a1 = Atom(name="CA", altloc="A", occupancy=0.3)
    a2 = Atom(name="CA", altloc="B", occupancy=0.7)
    a3 = Atom(name="CB")

    t = Topology()
    t.add_atom(a1)
    t.add_atom(a2)
    t.add_atom(a3)

    atoms = [a2, a3]  # ignore a1, lower occupancy.
    for idx, atom in enumerate(t):
        assert atom == atoms[idx]


def test_unpack_atoms():
    """Iterate over all atoms, incl altlocs, in the Topology."""

    a1 = Atom(name="CA", altloc="A")
    a2 = Atom(name="CA", altloc="B")
    a3 = Atom(name="CB", altloc="A")

    t = Topology()
    t.add_atom(a1)
    t.add_atom(a2)
    t.add_atom(a3)

    atoms = [a1, a2, a3]
    for idx, atom in enumerate(t.unpack_atoms()):
        assert atom == atoms[idx]


# Bonds
def test_set_bonds():
    """Set bondgraph in a Topology."""

    atoms = [Atom(name="CA"), Atom(name="CB"), Atom(name="CG"), Atom(name="CD")]

    t = Topology()
    for a in atoms:
        t.add_atom(a)

    with pytest.raises(TopologyError) as excinfo:
        t.bonds = []
    assert "Bond graph object must be a nx.Graph instance" in str(excinfo)

    with pytest.raises(TopologyError) as excinfo:
        t.bonds = nx.Graph()
    assert "Graph and Topology sizes do not match" in str(excinfo)

    t.bonds = nx.Graph([(0, 1), (1, 2), (2, 3)])  # OK!
