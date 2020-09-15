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

"""Unit tests for Atom class."""

import pytest

from interfacea.core.atom import Atom


def test_instantiate_error():
    """Throw error when creating empty Atom."""

    with pytest.raises(TypeError):
        _ = Atom()


def test_instantiate_success():
    """Create Atom (simple)."""

    a = Atom(name="CA")
    assert a.name == "CA"
    assert a.index is None


def test_instantiate_success_wmetadata():
    """Create Atom with metadata."""

    metadata = {"hetatm": False, "resname": "ARG", "chain": "A"}

    a = Atom(name="CA", **metadata)
    assert a.name == "CA"
    assert a.index is None
    assert a.hetatm is False
    assert a.resname == "ARG"
    assert a.chain == "A"


def test_atom_id():
    """Atom id property is set properly."""

    a = Atom(name="N")
    assert a.id == ("N", None, None, None)

    a = Atom(name="N", chain="A")
    assert a.id == ("N", "A", None, None)

    a = Atom(name="N", chain="A", resid=32)
    assert a.id == ("N", "A", 32, None)

    a = Atom(name="N", chain="A", resid=32, icode="A")
    assert a.id == ("N", "A", 32, "A")

    a = Atom(name="N", chain="A", resid=32, icode="A", altloc="A")
    b = Atom(name="N", chain="A", resid=32, icode="A", altloc="B")
    assert a.id == b.id


def test_atom_full_id():
    """Atom full_id property is set properly."""

    a = Atom(name="N")
    assert a.full_id == ("N", None, None, None, None)

    a = Atom(name="N", chain="A")
    assert a.full_id == ("N", "A", None, None, None)

    a = Atom(name="N", chain="A", resid=32)
    assert a.full_id == ("N", "A", 32, None, None)

    a = Atom(name="N", chain="A", resid=32, icode="A")
    assert a.full_id == ("N", "A", 32, "A", None)

    a = Atom(name="N", chain="A", resid=32, icode="A", altloc="D")
    assert a.full_id == ("N", "A", 32, "A", "D")


def test_atom_equality():
    """Compare two atoms."""

    a = Atom(name="N", chain="A", resid=32)
    b = Atom(name="N", chain="A", resid=32)
    c = Atom(name="C", chain="A", resid=32)

    assert a == b
    assert a is not b
    assert a != c
    assert not a == "string"


def test_atom_residue():
    """Print atom residue information (if available)."""

    a = Atom(name="N", chain="A", resname="GLY", resid=32, icode="A")
    assert a.residue == "Residue GLY32[A] of chain A"
    a.icode = None
    assert a.residue == "Residue GLY32 of chain A"
    a.chain = None
    assert a.residue == "Residue GLY32"
    a.chain = "A"
    a.resname = None
    assert a.residue == "Residue 32 of chain A"


def test_atom_copy():
    """Make copies of an Atom."""

    import copy

    a = Atom(name="N", chain="A", resid=32)
    shallow = copy.copy(a)

    assert a == shallow
    assert a is not shallow


def test_str_dunder():
    """Pretty printing works."""

    a = Atom(name="CA")
    assert a.__str__() == "<Atom CA index=None>"
