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
Unit tests for elements module.
"""

import interfacea.chemistry.elements as elements


def test_create_element():
    """Creates new element from base class"""

    e = elements.Element('C', 'Carbon', 6)
    assert e.symbol == 'C'
    assert e.name == 'Carbon'
    assert e.z == 6


def test_defined_elements():
    """Validate all known elements in ElementMapper"""

    for e in elements.ElementMapper.values():
        assert isinstance(e, elements.Element)


def test_known_element_lookup():
    """Retrieves element from ElementMapper"""

    e = elements.ElementMapper.get('C')
    assert e is elements.Carbon


def test_unknown_element_lookup():
    """Retrieves unknown element from ElementMapper"""

    e = elements.ElementMapper[0]
    assert e is elements.Unknown
