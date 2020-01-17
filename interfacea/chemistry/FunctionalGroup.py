#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018 João Pedro Rodrigues
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
Analysis of Biomolecular Interfaces.

Class(es) to define/match/find functional groups in molecules.
"""

import logging

import networkx as nx
import networkx.algorithms.isomorphism as iso

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


class FunctionalGroupError(Exception):
    """Custom Error Class"""
    pass


class FunctionalGroup(object):
    """Base class to represent chemical functional groups.

    Instances and subclasses can implement a `match` method
    and must define their own list of `elements` and `bonds`
    and `max_bonds`.

    Attributes:
        name (str): name of the functional group.
        charge (int): formal charge of the group.
        elements (`list(int)` or `list(tuple)`): elements contained in the group.
            Elements should be represented by their atomic number (1 - Hydrogen, 6 - Carbon, etc)
            A '0' (zero) is the wildcard, meaning any element matches.
            You can also pass a tuple of allowed elements, e.g. (1, 6) means match carbon or hydrogen.
        bonds (`list(tuple)`): bonds between elements of the functional group.
            Bond indexes refer to the positions of the elements in the `elements` set.
        max_bonds (`list(int)`): maximum number of bonds allowed for each element.
            Useful to restrict some searches (e.g. carboxylate vs esther). See examples below.
            Defaults to an arbitrary large number of bonds per atom (no filtering).

    Examples:

        >>> # Backbone (NH-CH-C-O)
        >>> #
        >>> #       H    H    O
        >>> #       |    |    |
        >>> #    -- N -- C -- C --
        >>> #            |
        >>> #            X
        >>> #
        >>> elements = [7, 1, 6, 1, 0, 6, 8]
        >>> bonds = [(0, 1), (0, 2), (2, 3), (2, 4), (2, 5), (5, 6)]

        >>> bb = fgs.FunctionalGroup(name='backbone',
                                     charge=0,
                                     elements=elements,
                                     bonds=bonds)


        >>> # Esther (X-O-O-X)
        >>> #
        >>> #         O
        >>> #        /
        >>> #    -- C -- O -- X
        >>> #
        >>> elements = [6, 8, 8, 0]
        >>> bonds = [(0, 1), (0, 2), (2, 3)]

        >>> esther = fgs.FunctionalGroup(name='esther',
                                         charge=0,
                                         elements=elements,
                                         bonds=bonds)

        >>> # Carboxylate (X-O-O)
        >>> #
        >>> #         O
        >>> #        /
        >>> #    -- C -- O
        >>> #
        >>> elements = [6, 8, 8]
        >>> bonds = [(0, 1), (0, 2)]
        >>> max_bonds = [3, 1, 1]

        >>> coo = fgs.FunctionalGroup(name='carboxylate',
                                      charge=0,
                                      elements=elements,
                                      bonds=bonds
                                      max_bonds=max_bonds)

    """

    __slots__ = ['name', 'charge', 'elements',
                 'bonds', 'max_bonds',
                 '_g', '_fixed_elem', '_fuzzy_elem']

    def __init__(self, name=None, charge=None, elements=None, bonds=None, max_bonds=None):

        assert name is not None, 'FunctionalGroup must have a name!'
        assert charge is not None, 'FunctionalGroup must have a charge!'
        assert elements is not None, 'FunctionalGroup must have atoms!'
        assert bonds is not None, 'FunctionalGroup must have bonds!'

        self.name = name
        self.elements = elements
        self.bonds = bonds
        self.charge = charge

        if max_bonds is None:
            max_bonds = [99 for e in elements]
        else:
            assert len(max_bonds) == len(elements), 'Items in max_bonds must match elements'
        self.max_bonds = max_bonds

        # Ensure all elements are bonded
        for idx, elem in enumerate(elements):
            in_bond = sum([1 for b in bonds if idx in b])
            if not in_bond:
                emsg = 'Atom #{} ({}) is not in any bond.'.format(idx, elem)
                raise FunctionalGroupError(emsg)

        # Ensure all bonds belong to existing elements
        for idx, bond in enumerate(bonds):
            a1, a2 = bond
            if a1 > len(elements) or a2 > len(elements):
                emsg = 'Bond #{} \'({}, {})\' includes an unknown atom.'
                raise FunctionalGroupError(emsg.format(idx, a1, a2))

        # Make element list a list of sets for efficient search
        _elements = []
        for elem in elements:
            if isinstance(elem, int):
                _elements.append({elem})
            elif isinstance(elem, (list, tuple)):
                _elements.append(set(elem))

        self.elements = _elements

        # Build list of fixed/fuzzy elements for matching
        # more efficiently when mathcing.
        self._fixed_elem = {next(iter(s)) for s in self.elements if len(s) == 1} - {0}
        self._fuzzy_elem = {e for s in self.elements for e in s if len(s) > 1} - {0}

        self._build_graph_representation()
        logging.debug('Successfully setup functional group \'{}\''.format(name))

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

        Does *not* match bonds across residues.

        Arguments:
            residue (:obj: `Residue`): OpenMM Residue object to scan for FG matches.

        Returns:
            All groups of atoms matching the FG in the `Residue`, as a list of lists of
            `Atom` objects.
        """

        matched_groups = []  # we can have more than one match!

        # Match atoms/elements first
        atomlist = list(residue.atoms())
        elemset = {a.element.atomic_number for a in atomlist}

        # Ensure all fixed elements are in the residue
        # and at least one fuzzy element is there too.
        no_fixed = not (elemset >= self._fixed_elem)  # subset
        no_fuzzy = bool(self._fuzzy_elem) and self._fuzzy_elem.isdisjoint(elemset)
        if no_fixed or no_fuzzy:
            return matched_groups

        # Match fg subgraph to residue graph
        def _node_match(n1, n2):
            """Return True if the residue atomic elements match the FG node elements.
            """
            n1_element = n1['element']
            n2_element = n2['element']
            return n1_element in n2_element or 0 in n2_element

        res_graph = residue._g
        matcher = iso.GraphMatcher(res_graph, self._g, node_match=_node_match)
        _mb = self.max_bonds  # make local variable
        if matcher.subgraph_is_isomorphic():

            observed = set()  # avoid degenerate matches
            for subgraph in matcher.subgraph_isomorphisms_iter():
                sorted_keys = tuple(sorted(subgraph.keys()))
                if sorted_keys in observed:
                    continue

                observed.add(sorted_keys)

                matched_atoms = []
                for r_idx in sorted_keys:
                    g_idx = subgraph[r_idx]
                    atom = atomlist[r_idx]

                    # Now check if atoms make only the bonds specified in the FG
                    num_bonds_in_res = len(residue.bonds_per_atom.get(atom))
                    if num_bonds_in_res > _mb[g_idx]:
                        logging.debug('Atom \'{}\' exceeds max_bonds ({}/{})'.format(atom.name,
                                                                                     num_bonds_in_res,
                                                                                     _mb[g_idx]))
                        break

                    matched_atoms.append(atom)
                else:
                    matched_groups.append(matched_atoms)

            logging.debug('{} matches functional group \'{}\' {} times'.format(residue,
                                                                               self.name,
                                                                               len(matched_groups)))
        return matched_groups

    def search(self, structure):
        """Utility function to search an entire structure matches of the functional group.

        Returns a dictionary with `Residue` objects as keys and the matched groups as values.
        """

        matches = {}
        for res in structure.topology.residues():
            matched_atoms = self.match(res)
            if matched_atoms:
                matches[res] = matched_atoms

        return matches
