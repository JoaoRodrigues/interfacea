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

"""Module containing classes to represent molecular topologies."""

import copy
import logging

import networkx as nx

from interfacea.core.atom import Atom, DisorderedAtom
from interfacea.exceptions import InterfaceaError


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class TopologyError(InterfaceaError):
    """Exception specific to submodule."""

    pass


class Topology:
    """Class to represent molecular topologies.

    The internal structure of a topology is a reduced version of what we
    see in mmCIF and PDB files, where molecular structures are divided in
    models, chains, residues, and finally atoms. In interfacea, topologies
    are bags of atoms. Disorder is handled at the atom level, with multiple
    occupancies bundling atoms in DisorderedAtoms. We do not handle point
    mutations (residues with the same sequence number, chain, and insertion
    code).

    Attributes
    ----------
        atoms : list
            ordered list of all atoms in the topology.
        bonds : networkX.Graph, optional
            bond graph of all atoms in the topology. Default is None.
        nlocs : int
            number of atoms in the topology, including alternate locations.
        natoms : int
            number of atoms in the topology, excluding alternate locations.
    """

    def __init__(self):
        """Create a new Topology."""

        logger.debug(f"Instantiated new Topology object: {id(self)}")

        self.natoms = 0
        self.nlocs = 0

        self.atoms = []
        self._bonds = None

        # Maps atom.id to topology indexes in atom list.
        self._atomdict = {}
        # Map loc index to atom index
        # DisorderedAtoms aggregate several locs in one atom.
        self._i2i = {}

    # Overloaded Dunders
    def __str__(self):
        """Return a string summary of the Topology."""

        return f"Topology with {self.natoms} atoms"

    def __copy__(self):
        """Return a shallow copy of the Topology."""

        new = Topology()
        for atom in self.atoms:
            new.add_atom(atom)

        if self.bonds is not None:
            new.bonds = self.bonds

        return new

    def __deepcopy__(self, memo):
        """Return a deep copy of the Topology."""

        new = Topology()
        for atom in self.atoms:
            new.add_atom(copy.deepcopy(atom, memo))

        if self.bonds is not None:
            new.bonds = copy.deepcopy(self.bonds, memo)

        return new

    def __getitem__(self, idx):
        """Retrieve an atom by its index in the Topology."""
        _idx = self._i2i.get(idx)
        if _idx is None:
            raise IndexError(f"Index '{idx}' not found in Topology")
        return self.atoms[_idx]

    def __setitem__(self, idx, atom):
        """Set an atom at a specific index in the Topology."""
        _idx = self._i2i.get(idx)
        if _idx is None:
            raise IndexError(f"Index '{idx}' not found in Topology")

        assert isinstance(atom, (Atom, DisorderedAtom))
        self.atoms[_idx] = atom

    def __iter__(self):
        """Iterate over children."""
        yield from self.atoms

    # Public Methods
    def add_atom(self, atom, disorder=True):
        """Add an Atom object to the Topology.

        Sets Atom.index based on Topology content.

        Parameters
        ----------
            atom : Atom
                individual atom to add to the topology, as an Atom instance.
            disorder : bool, optional
                if True, alternate conformations of the same atom are bundled
                in a DisorderedAtom object. If False, alternate locations are
                ignored and only the first instance is added to the topology.
        """

        id_in_topology = self._atomdict.get(atom.id)
        if id_in_topology is None:

            self.atoms.append(atom)

            # update counters and mapping
            self._atomdict[atom.id] = atom.index = self.nlocs  # set index
            self._i2i[self.nlocs] = self.natoms
            self.natoms += 1
            self.nlocs += 1

            logger.debug(f"Added atom '{atom}' at index {atom.index}")
        elif disorder:

            existing_atom = self[id_in_topology]
            if isinstance(existing_atom, Atom):
                da = DisorderedAtom()
                da.add(existing_atom)

                existing_atom = self[id_in_topology] = da
                logger.debug(f"New disordered atom at index {id_in_topology}")

            existing_atom.add(atom)

            # Update counter and mapping
            atom.index = self.nlocs  # set index
            self._i2i[self.nlocs] = self.natoms
            self.nlocs += 1
        else:
            og_atom = self[id_in_topology]
            logger.warning(f"Ignored alternate location of atom {og_atom}: {atom}")

    def unpack_atoms(self):
        """Return a generator of all atoms (including altlocs)."""

        for atom in self:
            if isinstance(atom, DisorderedAtom):
                for altloc in atom:
                    yield altloc
            else:
                yield atom

    @property
    def bonds(self):
        """Return a bond graph describing atom connectivity."""

        return self._bonds

    @bonds.setter
    def bonds(self, bg):
        """Set a bond graph describing atom connectivity.

        Arguments
        ---------
            bg : nx.Graph
                graph describing bonds between atoms in the structure. See
                interfacea.chemistry.bonds for details on the structure of the
                graph.
        """

        if not isinstance(bg, nx.Graph):
            raise TopologyError(
                f"Bond graph object must be a nx.Graph instance: {type(bg)}"
            )

        if len(bg) != self.natoms:
            raise TopologyError(
                f"Graph and Topology sizes do not match: {len(bg)} != {self.natoms}"
            )

        self._bonds = bg
        logger.debug(f"Set {bg} as topology bondgraph")
