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
Module containing es to define and search functional groups.
"""

import logging

import networkx as nx
import networkx.algorithms.components as nxcomp
import networkx.algorithms.isomorphism as nxiso

from interfacea.chemistry import elements as elem
from interfacea.exceptions import (
    FunctionalGroupError,
)

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    # Individual Groups
    'Carboxylate', 'Carboxyl', 'Guanidinium', 'Imidazole', 'Imidazolium',
    'Phosphate', 'SingleCoordinatedPhosphate', 'QuaternaryAmine', 'Sulfonium',
    'Sulfate', 'HydrogenSulfate', 'DivalentSulphur', 'AlkaneCarbon',
    'AlkeneCarbon', 'AlkyneCarbon', 'Phenyl', 'Indole', 'HBondDonor',
    # Collections
    'AllCharged', 'AllBases', 'AllAcids', 'AllHydrophobics', 'AllAromatics',
]

MAX_BONDS_PER_ATOM = 99


class FunctionalGroup(object):
    """Class to represent functional groups.

    A functional group is defined here as a collection of bonded atomic
    atomic elements. This simplistic representation is limited, but works well
    enough for most common groups in biomolecules.

    Arguments
    ----------
        name : str
            name of the functional group.

        elements : list
            Elements that make up the functional group. Use Unknown as a
            wildcard to match any element. You can also provide a tuples to
            represent multiple choices, e.g. (Oxygen, Nitrogen).

        bonds : list
            bonds between elements in the group, as tuples of numerical indices,
            e.g. (0, 1) is a bond between the first and second elements defined
            in the elements list above.

        max_bonds : list, optional
            maximum number of bonds that each element can make. We do this
            because we are unsure of bond orders and as such, cannot distinguish
            between C=O and C-O-X. Defaults to 99 bonds per element.

    Attributes
    ----------
        g : nx.Graph
            graph representation of the functional group, storing elements as
            nodes. The 'element' attribute is used by the `match` method to
            compare and find instances of g in a structure bond graph.
    """

    def __init__(self, name, elements, bonds, max_bonds=None):
        """Creates a new instance of the class"""

        self.g = None  # see _create_graph_representation

        self.name = name
        self.bonds = bonds

        if max_bonds is None:
            max_bonds = [MAX_BONDS_PER_ATOM for _ in elements]
        self.max_bonds = max_bonds

        assert len(elements) == len(max_bonds), \
            "Length of elements and max_bonds arguments must be equal."

        self.elements = [
            set((e,)) if isinstance(e, elem.Element) else set(e)
            for e in elements
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
        for idx, eset in enumerate(self.elements):
            g.add_node(idx, elem=eset)

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
        return (
            bg_node['elem'] in fg_node['elem'] or  # noqa:W504
            elem.Any in fg_node['elem']
        )

    def _exceeds_max_bonds(self, bondgraph, subgraph):
        """Returns True if any Atom exceeds the specified number of bonds."""

        for bg_idx, fg_idx in subgraph.items():
            if len(bondgraph.adj[bg_idx]) > self.max_bonds[fg_idx]:
                logging.debug(
                    f"Max bonds exceeded for atom #{bg_idx}"
                )
                return True
        return False

    def _is_missing_fg_elements(self, bondgraph):
        """Return True if any FG element is missing from the bondgraph.

        Arguments
        ---------
            bondgraph : nx.Graph
                bond graph calculated from a Structure object. Nodes must have
                an 'elem' attribute containing that atom's element as an Element
                object.
        """

        bondgraph_elementset = set(dict(bondgraph.nodes.data('elem')).values())
        bondgraph_elementset.update((elem.Any,))  # add wildcard

        for e in self.elements:
            if e.isdisjoint(bondgraph_elementset):
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

        logging.debug(f'Searching for instances of {self.name} in structure')

        matches = []

        # Get bond graph from structure
        bondgraph = structure.bonds

        if self._is_missing_fg_elements(bondgraph):
            return matches

        # Search for isomorphic subgraphs in structure
        m = nxiso.GraphMatcher(
            bondgraph,
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
                if self._exceeds_max_bonds(bondgraph, subgraph):
                    continue

                # Map to Atoms and add to matches
                matches.append(
                    [structure.atom(idx) for idx in subgraph]
                )

        logging.info(
            f"Found {len(matches)} matches of '{self.name}' in structure"
        )
        return matches


#
# Functional groups as singletons
#

# Charged
Carboxylate = FunctionalGroup(
    name='carboxylate',
    elements=[elem.Carbon, elem.Oxygen, elem.Oxygen],
    bonds=[
        (0, 1), (0, 2)
    ],
    max_bonds=[3, 1, 1]
)

Carboxyl = FunctionalGroup(
    name='carboxyl',
    elements=[
        elem.Carbon,
        elem.Oxygen,
        elem.Oxygen,
        elem.Hydrogen
    ],
    # the graph matching algorithm will
    # match either oxygens to the hydrogen
    # effectively having 'fuzzy' edges
    bonds=[
        (0, 1), (0, 2), (1, 3)
    ],
    max_bonds=[3, 1, 2, 1]

)

Guanidinium = FunctionalGroup(
    name='guanidinium',
    elements=[
        elem.Nitrogen,
        elem.Hydrogen,
        elem.Carbon,
        elem.Nitrogen,
        elem.Nitrogen,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Hydrogen
    ],
    bonds=[
        (0, 1), (0, 2), (2, 3), (2, 4),
        (3, 5), (3, 6), (4, 7), (4, 8)
    ],
    max_bonds=[3, 1, 3, 3, 3, 1, 1, 1, 1]
)

Imidazole = FunctionalGroup(
    name='imidazole',
    elements=[
        elem.Carbon,
        elem.Carbon,
        elem.Nitrogen,
        elem.Carbon,
        elem.Nitrogen
    ],
    bonds=[
        (0, 1), (1, 2), (2, 3), (3, 4),
        (4, 0)
    ],
    max_bonds=[3, 3, 3, 3, 3]
)

Imidazolium = FunctionalGroup(
    name='imidazolium',
    elements=[
        elem.Carbon,
        elem.Carbon,
        elem.Nitrogen,
        elem.Carbon,
        elem.Nitrogen,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Hydrogen
    ],
    bonds=[
        (0, 1), (1, 2), (1, 5), (2, 3),
        (2, 6), (3, 4), (3, 7), (4, 0),
        (4, 8)
    ],
    max_bonds=[3, 3, 3, 3, 3, 1, 1, 1, 1]
)

Phosphate = FunctionalGroup(
    name='phosphate',
    elements=[
        elem.Phosphorous,
        elem.Oxygen,
        elem.Oxygen,
        elem.Oxygen
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3)
    ],
    max_bonds=[4, 1, 1, 1]
)

SingleCoordinatedPhosphate = FunctionalGroup(
    name='phosphate-h',
    elements=[
        elem.Phosphorous,
        elem.Oxygen,
        elem.Oxygen,
        elem.Oxygen,
        elem.Any
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3), (1, 4)
    ],
    max_bonds=[4, 2, 1, 1, 99]
)

QuaternaryAmine = FunctionalGroup(
    name='quaternary_amine',
    elements=[
        elem.Nitrogen,
        elem.Any,
        elem.Any,
        elem.Any,
        elem.Any
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3), (0, 4)
    ]
)

Sulfonium = FunctionalGroup(
    name='sulfonium',
    elements=[
        elem.Sulfur,
        (elem.Carbon, elem.Hydrogen),
        (elem.Carbon, elem.Hydrogen),
        (elem.Carbon, elem.Hydrogen)
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3)
    ],
    max_bonds=[3, 4, 4, 4]
)

Sulfate = FunctionalGroup(
    name='sulfate',
    elements=[
        elem.Sulfur,
        elem.Oxygen,
        elem.Oxygen,
        elem.Oxygen
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3)
    ],
    max_bonds=[4, 1, 1, 1]
)

HydrogenSulfate = FunctionalGroup(
    name='sulfate-h',
    elements=[
        elem.Sulfur,
        elem.Oxygen,
        elem.Oxygen,
        elem.Oxygen,
        elem.Hydrogen
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3), (1, 4)
    ],
    max_bonds=[4, 2, 1, 1, 1]
)

# Hydrophobic
DivalentSulphur = FunctionalGroup(
    name='divalent-sulphur',
    elements=[
        elem.Sulfur,
        elem.Any,
        elem.Any
    ],
    bonds=[
        (0, 1), (0, 2)
    ],
    max_bonds=[2, 4, 4]
)


AlkaneCarbon = FunctionalGroup(
    name='alkane',
    elements=[
        elem.Carbon,
        (elem.Hydrogen, elem.Carbon),
        (elem.Hydrogen, elem.Carbon),
        (elem.Hydrogen, elem.Carbon),
        (elem.Hydrogen, elem.Carbon)
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3), (0, 4)
    ]
)


AlkeneCarbon = FunctionalGroup(
    name='alkene',
    elements=[
        elem.Carbon,
        (elem.Hydrogen, elem.Carbon),
        (elem.Hydrogen, elem.Carbon),
        (elem.Hydrogen, elem.Carbon)
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3)
    ],
    max_bonds=[3, 4, 4, 4]
)


AlkyneCarbon = FunctionalGroup(
    name='alkene',
    elements=[
        elem.Carbon,
        (elem.Hydrogen, elem.Carbon),
        (elem.Hydrogen, elem.Carbon)
    ],
    bonds=[
        (0, 1), (0, 2), (0, 3)
    ],
    max_bonds=[2, 4, 4]
)


Phenyl = FunctionalGroup(
    name='phenyl',
    elements=[
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Any,
        elem.Any,
        elem.Any,
        elem.Any,
        elem.Any,
        elem.Any
    ],
    bonds=[
        (0, 1), (1, 2), (2, 3), (3, 4),
        (4, 5), (5, 0), (0, 6), (1, 7),
        (2, 8), (3, 9), (4, 10), (5, 11)
    ]
)


Indole = FunctionalGroup(
    name='indole',
    elements=[
        elem.Nitrogen,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Carbon,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Any,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Hydrogen,
        elem.Hydrogen
    ],
    bonds=[
        (0, 1), (1, 2), (2, 3), (3, 4),
        (4, 5), (5, 6), (6, 7), (7, 8),
        (8, 0), (3, 8), (0, 9), (1, 10),
        (2, 11), (4, 12), (5, 13), (6, 14),
        (7, 15)
    ]
)


HBondDonor = FunctionalGroup(
    name='hbond-donor',
    elements=[
        (elem.Nitrogen, elem.Oxygen),
        elem.Hydrogen
    ],
    bonds=[
        (0, 1)
    ],
    max_bonds=[99, 1]
)


#
# Collections
#
AllAcids = [Carboxylate, Phosphate, SingleCoordinatedPhosphate, Sulfate]
AllBases = [Guanidinium, Imidazolium, QuaternaryAmine, Sulfonium]
AllCharged = AllBases + AllAcids

AllHydrophobics = [AlkaneCarbon, Phenyl, Indole, DivalentSulphur]
AllAromatics = [Phenyl, Indole]
