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

"""Module containing classes to build molecular structures from PDB files."""

import warnings

import interfacea.chemistry.elements as elements
from interfacea.exceptions import InterfaceaWarning


# Auxiliary Methods
def guess_atom_element(fullname):
    """Infer the atomic element from the atom name.

    Adapted from Biopython"s Atom.py code.

    Arguments
    ---------
        fullname : str
            the full 4-char string defined in the PDB file.

    Returns
    -------
        Element object or elements.Unknown if ambiguous/unsuccessful.

    """

    stripped_name = fullname.strip()

    if fullname[0].isalpha() and not fullname[2:].isdigit():
        putative_elem = stripped_name
    elif stripped_name[0].isdigit():  # e.g. 1HE2
        putative_elem = stripped_name[1]
    else:
        putative_elem = stripped_name[0]

    element = elements.ElementMapper.get(putative_elem.capitalize())
    if element is None:
        warnings.warn(
            f"Could not assign element for atom {repr(fullname)}",
            InterfaceaWarning,
        )
        element = elements.Unknown

    return element
