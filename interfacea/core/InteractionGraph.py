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

Module containing the InteractionGraph class for representing molecular
interfaces (and structures) as graphs of residues connected by interactions.
"""

import logging
import networkx as nx

from simtk.openmm.app import Residue as _Residue

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


# Error Classes
class InteractionGraphException(Exception):
    """Generic Exception for InteractionGraph objects"""
    pass


class DuplicatedNodeException(InteractionGraphException):
    """Thrown if trying to add an existing node to the InteractionGraph"""
    pass


class NodeNotFoundException(InteractionGraphException):
    """Thrown if trying to access a non-existing node in the InteractionGraph"""
    pass


class EdgeNotFoundException(InteractionGraphException):
    """Thrown if trying to access a non-existing edge in the InteractionGraph"""
    pass


# InteractionGraph Data Structure
class InteractionGraph(object):
    """Data structure to represent molecular interfaces as a queriable graph.

    Class implements a wrapper around networkx.MultiGraph.

    Args:
        name (str): descriptor to identify this InteractionGraph object.

    Attributes:
        num_residues (int): number of residues in the InteractionGraph.

    Example:

        >>> from interfacea.core.InteractionGraph import InteractionGraph
        >>> ...
        >>> g = InteractionGraph(name='IGraph_1brs')
        >>>
        >>> g.add_residue(resA)  # resA is an interfacea/OpenMM Residue
        >>> g.add_residue(resB)
        >>>
        >>> g.add_interaction(resA, resB, type='geometric', class='hbond')
        >>>
    """

    def __init__(self, name='generic graph'):
        """
        """
        self.name = name

        self._g = nx.MultiGraph()
        self._node_index = {}  # store residue->node index for fast(er) lookup

        self.num_residues = self._g.number_of_nodes()

    def __repr__(self):
        """String representation"""
        return f'InteractionGraph ({self.name}) ({self.num_residues})'

    def __contains__(self, item):
        """Returns True/False if a Residue/Pair of Residues is in the Graph"""

        if isinstance(item, _Residue):
            return self._has_residue(item)
        elif isinstance(item, tuple):
            try:
                i, j, a = item
            except ValueError:
                try:
                    i, j = item
                    a = None
                except ValueError as e:
                    emsg = f'Expected Residue or (Residue, Residue[, Attr])'
                    raise InteractionGraphException(emsg) from e
            return self._has_interaction(i, j, a)

    # Private Methods
    def _has_residue(self, residue):
        """Returns True if residue is in InteractionGraph"""
        return residue in self._node_index

    def _has_interaction(self, res_i, res_j, attributes=None):
        """Returns True if interaction is in InteractionGraph.

        If attributes dict is provided, returns True only if there is an
        interaction between res_i and res_j with matching attributes.
        """

        i = self._node_index.get(res_i)
        j = self._node_index.get(res_j)

        if i is None or j is None:  # one of the residue is missing
            return False

        _idx = (i, j)
        if attributes is None:  # simple check
            return _idx in self._g.edges(data=False)

        edges = self._g.get_edge_data(i, j)
        if edges is not None:
            for _, e_attr in edges.items():
                if e_attr == attributes:  # attribute check
                    return True
        return False

    def _get_residue_idx(self, residue):
        """Returns the index of the residue in the MultiGraph or dies trying"""
        idx = self._node_index.get(residue)
        if idx is None:
            emsg = f'Residue not found in InteractionGraph: {residue}'
            raise NodeNotFoundException(emsg)
        return idx

    # Public Methods
    def add_residue(self, residue, ignore_existing=True):
        """Adds a Residue to the InteractionGraph.

        Args:
            residue (:obj:`Residue`): residue object from OpenMM Topology.
            ignore_existing (bool): silently ignores if the residue already
                exists.
        Raises:
            DuplicatedNodeException: when trying to add an existing residue.
        """

        g = self._g

        if residue in self:
            node_idx = self._node_index.get(residue)
            if ignore_existing:
                return node_idx  # shhh

            emsg = f'{residue} is already in InteractionGraph: {node_idx}'
            raise DuplicatedNodeException(emsg)

        # new node
        node_idx = g.number_of_nodes() + 1
        g.add_node(node_idx, name=residue.name, id=residue.id,
                   icode=residue.insertionCode, chain=residue.chain.id)

        # Add to mapping
        self._node_index[residue] = node_idx
        logging.debug(f'Residue {residue} added to InteractionGraph: {node_idx}.')

        return node_idx

    def add_interaction(self, res_i, res_j, attributes=None):
        """Adds an interaction between two residues in the InteractionGraph.

        Args:
            res_i (:obj:`Residue`): residue object from OpenMM Topology.
            res_j (:obj:`Residue`): residue object from OpenMM Topology.
            attributes (dict): key/value pairs to add as attributes of the
                interaction. Default is None.

        Raises:
            DuplicatedEdgeException: when trying to add an existing interaction.
        """

        g = self._g
        attr = attributes

        idx_i = self._get_residue_idx(res_i)
        idx_j = self._get_residue_idx(res_j)

        # Multigraph might accumulate multiple edges with equivalent info
        has_edge = self._has_interaction(res_i, res_j, attributes)
        if has_edge:
            emsg = f'Skipping duplicated interaction {idx_i}-{idx_j}: {attr}'
            logging.info(emsg)

            # Get edge key
            edges = g.get_edge_data(idx_i, idx_j)
            for ekey, e_attr in edges.items():
                if e_attr == attr:  # attribute check
                    return ekey

        # New interaction
        ekey = g.add_edge(idx_i, idx_j)
        if attr is not None:
            g[idx_i][idx_j][ekey].update(attr)

        return ekey

    # Query Methods
    def get_residue_interactions(self, res_i, res_j, **kwargs):
        """Retrieves interaction(s) between residues i and j, optionally
        matching criteria.
        """

        g = self._g
        # Convert kwargs to set of tuples
        kw_set = {(k, v) for k, v in kwargs}

        # Multigraph might accumulate multiple edges with equivalent info
        has_edge = self._has_interaction(res_i, res_j)
        if not has_edge:
            emsg = f'Interaction not found between residues {res_i} and {res_j}'
            raise EdgeNotFoundException(emsg)

        # Get interactions between res_i and res_j
        idx_i = self._get_residue_idx(res_i)
        idx_j = self._get_residue_idx(res_j)

        result = []
        edges = g.get_edge_data(idx_i, idx_j)
        for ekey, e_attr in edges.items():
            if kw_set:
                if kw_set & {(k, v) for k, v in e_attr}:
                    result.append(e_attr)  # picky choice
            else:
                result.append(e_attr)

        return result
