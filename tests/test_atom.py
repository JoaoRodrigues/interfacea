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


def test_instantiate_error():
    """Throw error when creating empty Atom"""

    with pytest.raises(TypeError):
        _ = Atom()


def test_instantiate_success():
    """Create Atom (simple)"""

    a = Atom(name='CA')
    assert a.name == 'CA'
    assert a.serial is None
    assert a.parent is None

    with pytest.raises(AttributeError):
        _ = a.coords


def test_instantiate_success_wmetadata():
    """Create Atom with metadata"""

    metadata = {'hetatm': False, 'resname': 'ARG', 'chain': 'A'}

    a = Atom(name='CA', **metadata)
    assert a.name == 'CA'
    assert a.serial is None
    assert a.parent is None
    assert a.hetatm is False
    assert a.resname == 'ARG'
    assert a.chain == 'A'


def test_set_parent():
    """Sets Atom parent"""

    class Dummy:
        pass

    d = Dummy()

    a = Atom(name='N')
    a.parent = d
    assert a.parent is d


def test_del_parent():
    """Deletes Atom parent"""

    class Dummy:
        pass

    d = Dummy()

    a = Atom(name='N')
    a.parent = d
    del a.parent
    assert a.parent is None


def test_str_dunder():
    """Pretty printing works"""

    a = Atom(name='CA')
    assert a.__str__() == '<Atom name=CA serial=None>'