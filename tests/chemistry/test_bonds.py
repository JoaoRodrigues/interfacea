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
Unit tests for chemistry.bond submodul classes:
    - SimpleBondAnalyzer
"""

import itertools

from interfacea.io import read as io_read
from interfacea.chemistry.bonds import (
    SimpleBondAnalyzer
)

from . import datadir

defaultpdb = str(datadir / 'pdb' / 'triatom.pdb')


# SimpleBondAnalyzer
def test_SimpleBondAnalyzer():
    """Calculates bonds using SimpleBondAnalyzer"""

    s = io_read(defaultpdb)
    ba = SimpleBondAnalyzer(s)

    bonds = list(ba.run())
    assert len(bonds) == 1

    # intra-residue N-H bond only
    for aa, ab in itertools.permutations(s.atoms, 2):
        if aa.resid == ab.resid and aa.name != ab.name:
            assert aa in bonds[0] and ab in bonds[0]
        else:
            assert not (aa in bonds[0] and ab in bonds[0])


def test_SimpleBondAnalyzer_ignore():
    """SimpleBondAnalyzer ignores selected elements"""

    s = io_read(defaultpdb)
    ba = SimpleBondAnalyzer(s)
    ba.ignore = [1]  # ignore hydrogens

    bonds = list(ba.run())
    assert len(bonds) == 0
