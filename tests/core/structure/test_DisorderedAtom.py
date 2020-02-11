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
    - DisorderedAtom
"""

import pytest

from interfacea.core.structure import (
    Atom,
    DisorderedAtom,
)
import interfacea.exceptions as e

# Globals
aa = Atom('N', 0, altloc='A', occ=0.2)
ab = Atom('C', 1, altloc='B', occ=0.1)
ac = Atom('H', 2, altloc='C', occ=0.7)


# class DisorderedAtom
def test_create_DisorderedAtom():
    """Create DisorderedAtom from Atom objects."""

    da = DisorderedAtom()
    da.add(aa)
    assert da.nlocs == 1
    assert da.selected_child == aa
    assert da.name == 'N'

    da.add(ab)
    assert da.nlocs == 2
    assert da.selected_child == aa
    assert da.name == 'N'

    da.add(ac)
    assert da.nlocs == 3
    assert da.selected_child == ac
    assert da.name == 'H'


def test_create_DisorderedAtom_fromlist():
    """Populates DisorderedAtom using from_list method"""

    da = DisorderedAtom()
    da.from_list([aa, ab])

    assert da.nlocs == 2
    assert da.selected_child == aa
    assert da.name == 'N' and da.serial == 0


def test_empty_DisorderedAtom():
    """Raises error when accessing attributes of an empty DisorderedAtom."""

    da = DisorderedAtom()
    with pytest.raises(AttributeError):
        _ = da.name  # no children yet


def test_DisorderedAtom_DuplicateAltLocError():
    """Raises error when adding an existing altloc to DisorderedAtom"""

    da = DisorderedAtom()
    da.from_list([aa, ab])

    with pytest.raises(e.DuplicateAltLocError):
        da.add(aa)


def test_DisorderedAtom_select():
    """Changes selected_child on DisorderedAtom.select"""

    da = DisorderedAtom()
    da.from_list([aa, ab, ac])

    assert da.selected_child == ac
    assert da.name == 'H' and da.serial == 2

    da.select('B')

    assert da.selected_child == ab
    assert da.name == 'C' and da.serial == 1

    with pytest.raises(KeyError):
        da.select('Z')


def test_DisorderedAtom_setattr():
    """Forwards set-attribute calls to selected child"""

    da = DisorderedAtom()
    da.from_list([aa, ab, ac])

    da.name = 'O'
    assert da.name == 'O'
    assert da.selected_child.name == 'O'
