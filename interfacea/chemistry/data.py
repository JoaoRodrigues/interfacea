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
Module containing data tables with assorted information.
"""

# Covalent atom radii used to infer connectivity
# Maps atomic number to radii in Angstrom
# Values taken from CCDC and Roger Sayle's website:
# http://www.daylight.com/meetings/mug01/Sayle/m4xbondage.html
COVALENT_RADII = {
    1: 0.23,
    5: 0.83,
    6: 0.68,
    7: 0.68,
    8: 0.68,
    9: 0.64,
    14: 1.20,
    15: 1.05,
    16: 1.02,
    17: 0.99,
    33: 1.21,
    34: 1.22,
    35: 1.21,
    52: 1.47,
    53: 1.40
}
