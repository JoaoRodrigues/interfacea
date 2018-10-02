#!/usr/bin/env python

"""
Analysis of Biomolecular Interfaces.

Module containing definitions of functional groups
that can be matched to residues.
"""

from __future__ import print_function

import logging

import networkx as nx
import networkx.algorithms.isomorphism as iso

# Setup logger
# _private name to prevent collision/confusion with parent logger
_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


class FunctionalGroup(object):
    """Base class to represent chemical functional groups.

    Subclasses should implement a `match` method and have their
    own list of `elements` and `bonds`.

    e.g. defining an alcohol group

    oh_elements = [8, 1]
    oh_bonds = [(0, 1)]
    alcohol = FunctionalGroup(name='-oh', elements=oh_elements, bonds=oh_bonds)

    Attributes:
        name (str): name of the functional group.
        charge (int): formal charge of the group.
        elements (`list(int)`): atomic numbers of the elements contained in the group.
            A 0 (zero) acts as a wildcard, meaning any element.
        bonds (`list(tuple)`): bonds between elements of the functional group.
            Bond indexes refer to the positions of the elements in the `elements` set.
        terminal (`set(int)`): set of indexes of atoms with only one connection (terminal).
            Useful to distinguish between similar groups (e.g. esther vs carboxylate)
    """

    __slots__ = ['name', 'charge', 'elements', 'bonds', 'terminal']

    def __init__(self, name=None, charge=None, elements=None, bonds=None, terminal=None):

        assert name is not None, "FunctionalGroup must have a name!"
        assert charge is not None, "FunctionalGroup must have a charge!"
        assert elements is not None, "FunctionalGroup must have atoms!"
        assert bonds is not None, "FunctionalGroup must have bonds!"

        self.name = name
        self.elements = elements
        self.bonds = bonds
        self.charge = charge
        self.terminal = set([]) if terminal is None else set(terminal)

        self._element_set = set(elements)
        self._build_graph_representation()
        _logger.debug('Successfully setup functional group \'{}\''.format(name))

    def _build_graph_representation(self):
        """Builds a graph representation of the fragment.

        Creates a networkx Graph object with the atom of the
        functional group as nodes (element atomic number as
        an attribute) and bonds as edges.
        """

        g = nx.Graph()
        for idx, elem in enumerate(self.elements):
            g.add_node(idx, element=elem)

        for b in self.bonds:
            a1, a2 = b
            g.add_edge(a1, a2)

        self._g = g

    def match(self, residue):
        """Compares and returns matches between the FG and Residue graphs.

        Returns the atoms of the residue that match atoms of the functional
        group or an empty list otherwise.
        """

        # Match atoms/elements first
        atomlist = list(residue.atoms())
        elemlist = {a.element.atomic_number for a in atomlist}
        if not ((self._element_set - {0}) <= elemlist):
            return []  # not all atoms are in the Residue

        # Match fg subgraph to residue graph
        def _node_match(n1, n2):
            """Returns True if two nodes in the graph are equivalent.
            """
            # Edge matches, not node match.
            return (n1['element'] == n2['element'] or n2['element'] == 0)

        res_graph = residue._g
        matcher = iso.GraphMatcher(res_graph, self._g, node_match=_node_match)
        if matcher.subgraph_is_isomorphic():

            matched_atoms = []
            for r_idx, g_idx in matcher.mapping.items():
                atom = atomlist[r_idx]
                if g_idx in self.terminal and len(residue.bonds_per_atom.get(atom)) != 1:
                    return []
                matched_atoms.append(atom)

            return matched_atoms  # return the atoms matching the functional group

        return []


# Common functional groups
class Carboxylate(FunctionalGroup):
    """Carboxylate Functional Group.
    """

    def __init__(self):

        super().__init__(name='carboxylate',
                         charge=-1,
                         elements=[6, 8, 8],
                         bonds=[(0, 1), (0, 2)],
                         terminal={1, 2}
                         )


class Carboxyl(FunctionalGroup):
    """Carboxyl Functional Group.
    """

    def __init__(self):
        super().__init__(name='carboxyl',
                         charge=0,
                         elements=[6, 8, 8, 1],
                         # the graph matching algorithm will
                         # match both 8(1)-1 or 8(2)-1
                         # effectively having 'fuzzy' edges
                         bonds=[(0, 1), (0, 2), (1, 3)]
                         )


class Guanidinium(FunctionalGroup):
    """Guanidinium Functional Group.
    """

    def __init__(self):
        super().__init__(name='guanidinium',
                         charge=1,
                         elements=[6, 7, 7, 7, 1, 1, 1, 1],
                         bonds=[(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (3, 7)]
                         )


class Imidazole(FunctionalGroup):
    """Imidazole Functional Group.

    Without protons to avoid ambiguity between ND/NE protonation.
    User can check for position of proton later.
    """

    def __init__(self):
        super().__init__(name='imidazole',
                         charge=0,
                         elements=[6, 6, 7, 6, 7],
                         bonds=[(0, 1),
                                (1, 2),
                                (2, 3),
                                (3, 4),
                                (4, 0)]
                         )


class Imidazolium(FunctionalGroup):
    """Imidazolium Functional Group.
    """

    def __init__(self):
        super().__init__(name='imidazole',
                         charge=1,
                         elements=[6, 6, 7, 6, 7, 1, 1, 1, 1],
                         bonds=[(0, 1),
                                (1, 2), (1, 5),
                                (2, 3), (2, 6),
                                (3, 4), (3, 7),
                                (4, 0), (4, 8)]
                         )


class Phosphate(FunctionalGroup):
    """Phosphate Functional Group.
    """

    def __init__(self):
        super().__init__(name='phosphate',
                         charge=2,
                         elements=[15, 8, 8, 8],
                         bonds=[(0, 1), (0, 2), (0, 3)]
                         )


class HydrogenPhosphate(FunctionalGroup):
    """Hydrogen phosphate Functional Group.
    """

    def __init__(self):
        super().__init__(name='phosphate-h',
                         charge=1,
                         elements=[15, 8, 8, 8, 1],
                         bonds=[(0, 1), (0, 2), (0, 3), (1, 4)]
                         )


class QuaternaryAmine(FunctionalGroup):
    """Quaternary Amine Functional Group.
    """

    def __init__(self):
        super().__init__(name='quaternary_amine',
                         charge=1,
                         elements=[7, 1, 0, 0, 0],
                         bonds=[(0, 1), (0, 2), (0, 3), (0, 4)]
                         )


class Sulfate(FunctionalGroup):
    """Sulfate Functional Group.
    """

    def __init__(self):
        super().__init__(name='sulfate',
                         charge=1,
                         elements=[16, 8, 8, 8],
                         bonds=[(0, 1), (0, 2), (0, 3)]
                         )


class HydrogenSulfate(FunctionalGroup):
    """Hydrogen sulfate Functional Group.
    """

    def __init__(self):
        super().__init__(name='sulfate-h',
                         charge=0,
                         elements=[16, 8, 8, 8, 1],
                         bonds=[(0, 1), (0, 2), (0, 3), (1, 4)]
                         )
