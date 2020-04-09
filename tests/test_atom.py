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
Unit tests for Atom class.
"""

import pytest

from interfacea.core.atom import Atom


@pytest.fixture(scope='module')
def dummy():
    """Sets up a dummy class to act as a parent."""

    class Dummy:
        pass

    return Dummy()


def test_instantiate_error():
    """Throw error when creating empty Atom"""

    with pytest.raises(TypeError):
        _ = Atom()


def test_instantiate_success():
    """Create Atom (simple)"""

    a = Atom(name='CA')
    assert a.name == 'CA'
    assert a.serial == 0
    assert a.parent is None

    with pytest.raises(AttributeError):
        _ = a.coords


def test_instantiate_success_wmetadata():
    """Create Atom with metadata"""

    metadata = {'hetatm': False, 'resname': 'ARG', 'chain': 'A'}

    a = Atom(name='CA', **metadata)
    assert a.name == 'CA'
    assert a.serial == 0
    assert a.parent is None
    assert a.hetatm is False
    assert a.resname == 'ARG'
    assert a.chain == 'A'


def test_set_parent(dummy):
    """Sets Atom parent"""

    a = Atom(name='N')
    a.parent = dummy
    assert a.parent is dummy


def test_del_parent(dummy):
    """Deletes Atom parent"""

    a = Atom(name='N')
    a.parent = dummy
    del a.parent
    assert a.parent is None


def test_str_dunder():
    """Pretty printing works"""

    a = Atom(name='CA')
    assert a.__str__() == '<Atom name=CA serial=0>'
