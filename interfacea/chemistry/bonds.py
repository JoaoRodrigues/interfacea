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
Module containing classes to detect atom connectivity/bonds.
"""

import abc
import logging

import networkx as nx

from interfacea.chemistry.data import COVALENT_RADII

logging.getLogger(__name__).addHandler(logging.NullHandler())


class BaseBondDetector(abc.ABC):
    """Base class for bond detectors.

    Subclasses must implement a run() method that takes a Structure object
    as an argument and returns a networkx Graph defining bonds between all atoms
    in the structure.
    """

    @abc.abstractmethod
    def run(structure):
        """Analyzes a structure to infer atom connectivity.

        Arguments
        ---------
            structure : Structure
                atomic structure with coordinate and metadata.

        Returns
        -------
            a networkx Graph object representing the bond graph.
        """
        pass


class SimpleBondDetector(BaseBondDetector):
    """Infers atom connectivity from xyz coordinates. Very simplistic."""

    def __init__(self, tolerance=0.45, ignorelist=None):
        """
        Arguments
        ---------
            tolerance : float, optional
                fudge factor (in Angstrom) to account for variations in bond
                lengths. A lower value means a more strict criterion to consider
                two atoms bonded, so less bonds inferred. Default is 0.45.
            ignorelist : list, optional
                atomic numbers to ignore when accounting for bonding. By default
                contains all metallic elements.

        Raises
        ------
            ValueError
                if the tolerance value is zero or a negative number.
        """

        if tolerance <= 0:
            raise ValueError(f'Tolerance factor must be > 0')

        if ignorelist is None:
            ignorelist = [
                3, 4, 12, 13, 18, 23, 24, 25, 27, 29, 30, 31, 33, 36, 37, 38, 39,
                42, 47, 48, 49, 51, 52, 54, 55, 56, 57, 58, 62, 63, 64, 65, 67,
                70, 71, 74, 75, 76, 77, 78, 79, 80, 81, 82, 92,
            ]

        self.ignoreset = set(ignorelist)
        self.tolerance = tolerance

    def run(self, structure):
        """
        Arguments
        ---------
            structure : Structure
                atomic structure with coordinate and metadata.

        Returns
        -------
            a networkx Graph object representing the bond graph.
        """

        covradii = COVALENT_RADII

        # Initialize bond graph with atoms as nodes
        # and atomic numbers as attributes
        bondgraph = nx.Graph()
        for atom in structure:
            bondgraph.add_node(atom.serial, Z=atom.element.atomic_number)

        max_r = (max(covradii.values()) ** 2) + self.tolerance

        seen = set()  # store atoms we've visited not to repeat.
        for atom in structure:
            if atom.element.atomic_number in self.ignoreset:
                logging.debug(f"Skipping bonds for {atom}: ignored element")
                continue

            atom_r = covradii.get(
                atom.element.atomic_number,
                2.0
            )

            nlist = structure.neighbors(atom, radius=max_r)
            for neighbor, d in nlist:
                if (
                    neighbor in seen or  # noqa: W504
                    neighbor.element.atomic_number in self.ignoreset
                ):
                    continue

                neighbor_r = covradii.get(
                    neighbor.element.atomic_number,
                    2.0
                )
                if d <= (atom_r + neighbor_r + self.tolerance):
                    bondgraph.add_edge(atom.serial, neighbor.serial)

            seen.add(atom)

        return bondgraph
