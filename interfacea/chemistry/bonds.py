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

    Subclasses must implement a run() method that returns a networkx Graph
    defining bonds between (all) atoms in the structure. Nodes (atoms) must have
    at least an 'elem' attribute with the corresponding Element object from each
    atom. Edges are not required to have any attribute.
    """

    @abc.abstractmethod
    def run(self):
        """Analyzes a structure to infer atom connectivity.

        Returns
        -------
            nx.Graph object representing the bond graph, where nodes represent
            atoms and edges bonds between them.
        """
        pass


class SimpleBondDetector(BaseBondDetector):
    """Infers atom connectivity from xyz coordinates. Very simplistic."""

    def __init__(self, structure, tolerance=0.45, ignorelist=None):
        """
        Arguments
        ---------
            structure : Structure
                atomic structure with coordinates and metadata.

            tolerance : float, optional
                fudge factor (in Angstrom) to account for variations in bond
                lengths. A lower value means a more strict criterion to consider
                two atoms bonded, so less bonds inferred. Default is 0.45.
            ignorelist : list, optional
                elements to ignore when accounting for bonding, as a list of
                Element objects. By default contains all metallic elements.
        """

        IGNORED_ELEMENTS_Z = {
            3, 4, 12, 13, 18, 23, 24, 25, 27, 29, 30, 31, 33, 36, 37, 38, 39,
            42, 47, 48, 49, 51, 52, 54, 55, 56, 57, 58, 62, 63, 64, 65, 67,
            70, 71, 74, 75, 76, 77, 78, 79, 80, 81, 82, 92,
        }

        if tolerance <= 0:
            raise ValueError(f'Tolerance factor must be > 0')

        if ignorelist is None:
            ignoreset = IGNORED_ELEMENTS_Z
        else:
            try:
                ignoreset = {e.z for e in ignorelist}
            except TypeError as err:
                raise TypeError(
                    f"Could not convert ignore list to a set: {err}"
                )
            except AttributeError as err:
                raise AttributeError(
                    f"Items in ignore list must be Element instances: {err}"
                )

        self.structure = structure
        self.ignoreset = ignoreset
        self.tolerance = tolerance

        self.init_graph()

    def init_graph(self):
        """Initializes the bond graph with atoms from the Structure"""

        self.bondgraph = bondgraph = nx.Graph()

        for atom in self.structure:
            if atom.element.z in self.ignoreset:
                continue

            bondgraph.add_node(
                atom.serial,
                elem=atom.element
            )

        logging.debug(
            f"Initialized bond graph with {len(self.bondgraph)} nodes."
        )

    def run(self):
        """
        Returns
        -------
            nx.Graph object representing the bond graph, where atoms are the
            nodes and bonds between them the edges. Nodes have an 'elem'
            attribute corresponding to their respective atom.Element object.
        """

        covradii = COVALENT_RADII
        max_r = (max(covradii.values()) ** 2) + self.tolerance

        visited = set()  # avoid i-j, j-i dupes
        for idx in self.bondgraph:

            atom = self.structure.atom(idx)

            atom_r = covradii.get(atom.element.z, 2.0)

            visited.add(atom.serial)  # avoid i-i

            neighbors = self.structure.neighbors(atom, radius=max_r)
            for n, dist in neighbors:
                if n in visited or n.element in self.ignoreset:
                    continue

                neighbor_r = covradii.get(n.element.z, 2.0)
                if dist <= (atom_r + neighbor_r + self.tolerance):
                    self.bondgraph.add_edge(atom.serial, n.serial)

        return self.bondgraph
