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

from interfacea.exceptions import InterfaceaWarning
from interfacea.io import read
from interfacea.io.helpers import guess_atom_element

import pytest


# 1ctf: monomer, single-chain, no altlocs or icodes.
def test_reading_1ctf_as_str(datadir):
    """Read 1ctf.pdb as string."""

    read(str(datadir / "1ctf.pdb"))


def test_reading_1ctf_as_path(datadir):
    """Read 1ctf.pdb as path."""

    read(datadir / "1ctf.pdb")


def test_reading_1ctf_gz_as_path(datadir):
    """Read 1ctf.pdb.gz as path."""

    read(datadir / "1ctf.pdb.gz")


def test_reading_1ctf_with_topology(datadir):
    """Triggers warning when parsing 1ctf.pdb with topology argument."""

    with pytest.warns(InterfaceaWarning, match="Ignoring provided topology"):
        read(datadir / "1ctf.pdb", topology=str(datadir / "1ctf.pdb"))


# Helpers
def test_guess_element_from_atom_name():
    """Guess element from atom name."""

    answers = {
        " CA ": "C",
        "CA  ": "X",  # until we add calcium, unknown maps to X
        " N  ": "N",
        "HH21": "H",
        "1HE1": "H",
    }

    for name, element in answers.items():
        if element == "X":
            with pytest.warns(InterfaceaWarning, match="Could not assign"):
                e = guess_atom_element(name)
        else:
            e = guess_atom_element(name)
        assert e.symbol == element
