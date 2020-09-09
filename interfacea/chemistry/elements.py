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

"""Module containing data and classes to handle atomic elements."""

import collections


__all__ = [
    "Hydrogen",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Phosphorous",
    "Sulfur",
    "Unknown",
]


Element = collections.namedtuple("Element", ["symbol", "name", "z"])

# Singletons
Unknown = Element("X", "unknown", None)
Any = Element("*", "any", None)

Hydrogen = Element("H", "hydrogen", 1)
Carbon = Element("C", "carbon", 6)
Nitrogen = Element("N", "nitrogen", 7)
Oxygen = Element("O", "oxygen", 8)
Phosphorous = Element("P", "phosphorous", 15)
Sulfur = Element("S", "sulfur", 16)

# Mapping
ElementMapper = {}

_knowns = {
    "H": Hydrogen,
    "C": Carbon,
    "N": Nitrogen,
    "O": Oxygen,
    "P": Phosphorous,
    "S": Sulfur,
}

ElementMapper.update(_knowns)
