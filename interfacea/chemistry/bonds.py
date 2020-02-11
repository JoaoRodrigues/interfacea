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
Module containing bonding information and analyzers.
"""

import logging

from interfacea.core.structure import Structure
from interfacea.chemistry.data import COVALENT_RADII

logging.getLogger(__name__).addHandler(logging.NullHandler())


class SimpleBondAnalyzer(object):
    """Class to infer atom connectivity from cartesian coordinates.

    Args:
        structure (Structure): atomic structure with coordinate and metadata.

    Attributes:
        tolerance (float): factor to account for variations in bond lengths.
            Set to 0.45 by default, lower means a more strict criterion to
            consider two atoms bonded.
        ignore (list): atomic numbers to ignore when accounting for bonding.
            By default, contains all metallic elements.

    Returns:
        List of 2-item tuples with Atoms involved in a bond.
    """

    def __init__(self, structure):

        self.tolerance = 0.45
        self.ignore = [
            3, 4, 12, 13, 18, 23, 24, 25, 27, 29, 30, 31, 33, 36, 37, 38, 39,
            42, 47, 48, 49, 51, 52, 54, 55, 56, 57, 58, 62, 63, 64, 65, 67,
            70, 71, 74, 75, 76, 77, 78, 79, 80, 81, 82, 92,
        ]

        self._load_structure(structure)  # sets self.structure

    def _load_structure(self, structure):
        """Validates the input structure object

        Raises:
            TypeError: if input structure is not a Structure type.
        """

        if not isinstance(structure, Structure):
            emsg = f"Object is not a Structure (sub)class: {type(structure)}"
            raise TypeError(emsg)

        self.structure = structure

    def run(self):
        """Runs the analyzer."""

        seen = set()  # store bonds we've seen not to repeat.

        ignore_set = set(self.ignore)

        # Load kdtree
        if self.structure._kdtree:
            self.structure.generate_kdtree()

        max_r = (max(COVALENT_RADII.values()) ** 2) + self.tolerance
        get_r = lambda a: COVALENT_RADII.get(a.element.atomic_number, 2.0)

        for atom in self.structure:
            if atom.element.atomic_number in ignore_set:
                logging.debug(f"Skipping bonds for {atom}: ignored element")
                continue

            atom_r = get_r(atom)

            nlist = atom.neighbors(radius=max_r)
            for neighbor, d in nlist:
                if (
                    neighbor in seen or  # noqa: W504
                    neighbor.element.atomic_number in ignore_set
                ):
                    continue

                neighbor_r = get_r(neighbor)
                if d <= (atom_r + neighbor_r + self.tolerance):
                    yield (atom, neighbor)

            seen.add(atom)
