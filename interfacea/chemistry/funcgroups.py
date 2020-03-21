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
Module containing classes to define and search functional groups.
"""

import logging

import networkx as nx
import networkx.algorithms.components as nxcomp
import networkx.algorithms.isomorphism as nxiso

from interfacea.exceptions import (
    FunctionalGroupError,
)

logging.getLogger(__name__).addHandler(logging.NullHandler())


class BaseFunctionalGroup(object):
    """Base class to represent functional group.

    A functional group is defined here as a collection of bonded atomic
    atomic elements. This simplistic representation is limited, but works well
    enough for most common groups in biomolecules.

    Arguments
    ----------
        name : str
            name of the functional group.

        elements : list
            atomic numbers of the elements that make up the functional group.
            Use 0 as a wildcard to match any element. You can also provide a
            tuple of numbers to represent multiple choices, e.g. (8, 7) matches
            oxygen or nitrogen.

        bonds : list
            bonds between elements in the group, as tuples of numerical indices,
            e.g. (0, 1) is a bond between element 0 and 1 as defined in the
            elements list.

        max_bonds : list, optional
            maximum number of bonds that each element can make. We do this
            because we are unsure of bond orders and as such, cannot distinguish
            between C=O and C-O-X.

    Attributes
    ----------
        g : nx.Graph
            graph representation of the functional group, storing atoms as nodes
            with a 'Z' attribute containing a set of allowed atomic numbers.
            This attribute is used by the `match` method to compare and find
            instances of g in a structure bond graph.

        elementset : set
            list of sets containing element composition of the functional group.
    """

    def __init__(self, name, elements, bonds, max_bonds=None):
        """Creates a new instance of the class"""

        self.g = None  # see _create_graph_representation

        self.name = name
        self.bonds = bonds

        if max_bonds is None:
            max_bonds = [99 for _ in elements]
        self.max_bonds = max_bonds

        assert len(elements) == len(max_bonds), \
            "Length of elements and max_bonds arguments must be equal."

        self.elementset = [
            set((e,)) if isinstance(e, int) else set(e) for e in elements
        ]

        try:
            self._create_graph_representation()
        except Exception as err:
            raise FunctionalGroupError(
                f"Error when creating functional group '{name}': {err}"
            ) from err

        logging.debug(f"Successfully setup new instance of '{self.name}'")

    # Graph Creation
    def _create_graph_representation(self):
        """Converts the elements/bonds into a valid graph representation."""

        g = nx.Graph()
        for idx, eset in enumerate(self.elementset):
            g.add_node(idx, Z=eset)

        for e1, e2 in self.bonds:
            g.add_edge(e1, e2)

        if not nxcomp.is_connected(g):
            raise ValueError(
                f"Functional group '{self.name}' contains isolated elements"
            )

        self.g = g

    # Matching Methods
    def _node_matcher(self, bg_node, fg_node):
        """Returns True if a fg node matches a structure node."""
        return bg_node['Z'] in fg_node['Z'] or 0 in fg_node['Z']

    def _map_subgraph_to_structure(self, subgraph):
        """Returns a dictionary mapping Atom objects to FG nodes"""

        return [
            self.structure.atom(idx) for idx in subgraph
        ]

    def _exceeds_max_bonds(self, subgraph):
        """Checks if any of the matched Atoms exceeed its max_bonds parameter"""

        for bg_idx, fg_idx in subgraph.items():
            if len(self.bondgraph.adj[bg_idx]) > self.max_bonds[fg_idx]:
                logging.debug(
                    f"Max bonds exceeded for atom #{bg_idx}"
                )
                return True
        return False

    def _is_missing_fg_elements(self):
        """Return True if any FG element is missing from the structure."""

        struct_elements = {
            z for _, z in self.bondgraph.nodes.data('Z')
        }
        struct_elements.update((0,))  # add wildcard

        for eset in self.elementset:
            if eset.isdisjoint(struct_elements):
                logging.debug(
                    f"Elements of '{self.name}' missing in structure"
                )
                return True
        return False

    def match(self, structure):
        """Returns matches of the functional group in a Structure.

        Arguments
        ---------
            structure : Structure
                atomic structure to search for functional group.

        Returns
        -------
            list of lists, each containing a group of atoms that matches the
            functional group.
        """

        matches = []

        # Get bond graph from structure
        self.structure = structure
        self.bondgraph = structure.bonds

        if self._is_missing_fg_elements():
            return matches

        # Search for isomorphic subgraphs in structure
        m = nxiso.GraphMatcher(
            self.bondgraph,
            self.g,
            node_match=self._node_matcher
        )
        if m.subgraph_is_isomorphic():

            seen = set()  # avoid degenerate matches
            for subgraph in m.subgraph_isomorphisms_iter():

                indexes = tuple(sorted(subgraph))
                if indexes in seen:
                    continue
                seen.add(indexes)

                # Check connectivity for max_bonds violations
                if self._exceeds_max_bonds(subgraph):
                    continue

                # Map to Atoms
                subgraph_atoms = self._map_subgraph_to_structure(
                    subgraph
                )

                matches.append(subgraph_atoms)

        logging.info(
            f"Found {len(matches)} matches of '{self.name}' in structure"
        )
        return matches
