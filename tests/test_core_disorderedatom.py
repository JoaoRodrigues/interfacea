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

from unittest.mock import Mock

import pytest

from interfacea.core.atom import (
    DisorderedAtom,
    DisorderedAtomError
)


def make_mock_atom(**kwargs):
    """Make a mock Atom object with kwargs params.

    Example
    -------
        >>> a = make_mock_atom(name='CA')
        >>> assert a.name == 'CA'
        >>> with pytest.raises(AttributeError):
        >>>     _ = a.occ
    """
    m = Mock(spec=[])  # raises AttributeError on what we do not define below
    m.configure_mock(**kwargs)
    return m


def test_instantiate_success():
    """Intantiate the DisorderedAtom class"""

    da = DisorderedAtom()
    assert hasattr(da, 'children')
    assert hasattr(da, 'representative')


def test_str_dunder():
    """Pretty-printing of an empty DisorderedAtom"""

    da = DisorderedAtom()
    assert da.__str__() == '<Empty DisorderedAtom>'


def test_str_dunder_nonempty():
    """Pretty-printing of a DisorderedAtom"""

    da = DisorderedAtom()
    # override __setattr__
    da.__dict__['name'] = 'Hello'
    da.__dict__['serial'] = 0
    assert da.__str__() == "<DisorderedAtom name=Hello serial=0 nlocs=0>"


def test_add_atom():
    """Add Atoms to a DisorderedAtom"""

    a1 = make_mock_atom(name='CA', altloc='A', occ=0.7)
    a2 = make_mock_atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    assert da.occ == 0.7
    assert da.nlocs == 1
    da.add(a2)
    assert da.occ == 0.7
    assert da.nlocs == 2
    assert da.representative == a1


# see: https://docs.pytest.org/en/latest/logging.html#caplog-fixture
def test_add_atom_woutocc(caplog):
    """Add Atoms to a DisorderedAtom (no occ)"""

    a1 = make_mock_atom(name='CA', altloc='A')
    a2 = make_mock_atom(name='CA', altloc='B')

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

    a1 = make_mock_atom(name='CA', altloc='A', occ=0.3)  # swapped
    a2 = make_mock_atom(name='CA', altloc='B', occ=0.7)

    da = DisorderedAtom()
    da.from_list([a1, a2])
    assert da.occ == 0.7
    assert da.nlocs == 2


def test_autoset_altloc(caplog):
    """Auto-sets altloc (and issues warning) when adding atoms"""

    a1 = make_mock_atom(name='CA', occ=0.3)  # set occ not to trigger additional warning
    a2 = make_mock_atom(name='CA', occ=0.7)

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

    a1 = make_mock_atom(name='CA', altloc='A')

    da = DisorderedAtom()
    da.add(a1)
    with pytest.raises(DisorderedAtomError):
        da.add(a1)


def test_iterator():
    """Iterate over list of child atoms"""

    a1 = make_mock_atom(name='CA', altloc='A', occ=0.7)
    a2 = make_mock_atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)

    assert list(da) == list(da.children.values())


def test_set_child_attribute():
    """Set attribute of selected child"""

    a1 = make_mock_atom(name='CA', altloc='A', occ=0.7)
    a2 = make_mock_atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)

    da.name = 'N'  # stupid example
    assert da.name == da.representative.name == 'N'


def test_get_child_attribute_error():
    """Throws error when attribute of child is not found"""

    a1 = make_mock_atom(name='CA', altloc='A')
    da = DisorderedAtom()
    da.add(a1)

    with pytest.raises(AttributeError):
        _ = da.occ


def test_get_child_attribute_error_2():
    """Throws error when accessing attribute without children"""

    da = DisorderedAtom()

    with pytest.raises(DisorderedAtomError):
        _ = da.occ


def test_select_child():
    """Select child manually"""

    a1 = make_mock_atom(name='CA', altloc='A', occ=0.7)
    a2 = make_mock_atom(name='CA', altloc='B', occ=0.3)

    da = DisorderedAtom()
    da.add(a1)
    da.add(a2)

    assert da.representative == a1
    da.select("B")
    assert da.representative == a2


def test_select_child_error():
    """Throw error when selecting unknown altloc"""

    a1 = make_mock_atom(name='CA', altloc='B', occ=0.7)

    da = DisorderedAtom()
    da.add(a1)

    with pytest.raises(DisorderedAtomError):
        da.select("A")


def test_select_child_error2():
    """Throw error when selecting from empty DisorderedAtom"""

    da = DisorderedAtom()
    with pytest.raises(DisorderedAtomError):
        da.select("A")
