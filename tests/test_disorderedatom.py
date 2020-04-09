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
Unit tests for DisorderedAtom class.
"""

import pytest

from interfacea.core.atom import Atom, DisorderedAtom, DisorderedAtomError


def test_instantiate_success():
    """Intantiate the DisorderedAtom class"""

    da = DisorderedAtom()
    assert hasattr(da, 'children')
    assert hasattr(da, 'selected_child')


def test_str_dunder():
    """print(atom) does not raise an error"""

    da = DisorderedAtom()
    assert da.__str__() == '<Empty DisorderedAtom>'


def test_add_atom():
    """Add Atoms to a DisorderedAtom"""

    a1 = Atom(name='CA', altloc='A', occ=0.7)
    a2 = Atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    assert da.occ == 0.7
    assert da.nlocs == 1
    da.add(a2)
    assert da.occ == 0.7
    assert da.nlocs == 2
    assert da.selected_child == a1


# see: https://docs.pytest.org/en/latest/logging.html#caplog-fixture
def test_add_atom_woutocc(caplog):
    """Add Atoms to a DisorderedAtom (no occ)"""

    a1 = Atom(name='CA', altloc='A')
    a2 = Atom(name='CA', altloc='B')

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)  # should trigger log.warning

    assert len(caplog.records) == 1
    for record in caplog.records:
        assert record.levelname == "WARNING"

    assert da.altloc == 'B'  # last one added is selected
    assert da.nlocs == 2


def test_from_list():
    """Add a list of Atoms to a DisorderedAtom"""

    a1 = Atom(name='CA', altloc='A', occ=0.3)  # swapped
    a2 = Atom(name='CA', altloc='B', occ=0.7)

    da = DisorderedAtom()
    da.from_list([a1, a2])
    assert da.occ == 0.7
    assert da.nlocs == 2


def test_autoset_altloc(caplog):
    """Auto-sets altloc (and issues warning) when adding atoms"""

    a1 = Atom(name='CA', occ=0.3)  # set occ not to trigger additional warning
    a2 = Atom(name='CA', occ=0.7)

    da = DisorderedAtom()
    da.add(a1)
    assert da.altloc == "0"
    da.add(a2)
    assert da.altloc == "1"

    assert len(caplog.records) == 2
    for record in caplog.records:
        assert record.levelname == "WARNING"


def test_duplicated_altloc_error():
    """Throw error when adding duplicated altloc"""

    a1 = Atom(name='CA', altloc='A')

    da = DisorderedAtom()
    da.add(a1)
    with pytest.raises(DisorderedAtomError):
        da.add(a1)


def test_iterator():
    """Iterate over list of child atoms"""

    a1 = Atom(name='CA', altloc='A', occ=0.7)
    a2 = Atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)

    assert list(da) == list(da.children.values())


def test_set_child_attribute():
    """Set attribute of selected child"""

    a1 = Atom(name='CA', altloc='A', occ=0.7)
    a2 = Atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)

    da.name = 'N'  # stupid example
    assert da.name == da.selected_child.name == 'N'


def test_get_child_attribute_error():
    """Throws error when attribute of child is not found"""

    a1 = Atom(name='CA', altloc='A')
    da = DisorderedAtom()
    da.add(a1)

    with pytest.raises(AttributeError):
        _ = da.occ


def test_select_child():
    """Select child manually"""

    a1 = Atom(name='CA', altloc='A', occ=0.7)
    a2 = Atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)

    assert da.selected_child == a1
    da.select("B")
    assert da.selected_child == a2


def test_select_child_error():
    """Throw error when selecting unknown altloc"""

    da = DisorderedAtom()
    with pytest.raises(DisorderedAtomError):
        da.select("A")
