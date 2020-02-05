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
Contains AtomRecord class to store information read from coordinate files.
"""

from dataclasses import dataclass

import interfacea.chemistry.elements as e


# Data class to store atomic data _before_ building a structure.
@dataclass
class AtomRecord:
    """Stores data on an ATOM/HETATM line."""

    serial: int
    model: int
    rectype: str
    name: str
    altloc: str
    resname: str
    chain: str
    resid: int
    icode: str
    x: float
    y: float
    z: float
    occ: float
    b: float
    segid: str
    element: e.Element = e.unk  # default to unknown
