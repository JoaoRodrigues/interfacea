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
Module containing classes to represent 3D molecular objects.

Atomic structures in interfacea are represented internally by the Structure
class, which stores both coordinate data and metadata. The metadata, e.g. atom
and residue names, is stored in individual Atom objects.

Structures are meant to be manipulated but _not created_ by users.
"""

import logging

import networkx as nx
import numpy as np

from interfacea.core.atom import Atom, DisorderedAtom
from interfacea.chemistry.bonds import SimpleBondDetector
from interfacea.exceptions import InterfaceaError
from interfacea.spatial import kdtrees

logging.getLogger(__name__).addHandler(logging.NullHandler())


class StructureError(InterfaceaError):
    pass


class Structure(object):
    """Represents 3D molecules as collections of atoms.

    Attributes
    ----------
        name : str
            string that identifies the structure.
        atoms : list
            ordered list of all atoms in the structure.
        bonds : networkX.Graph
            bond graph of all atoms in the structure.
        coords : np.ndarray
            coordinate array of the selected model.
        raw_coords : np.ndarray
            raw coordinate array including data for all models.
        natoms : int
            number of atoms in the structure.
        nmodels : int
            number of models in the structure.
    """

    def __init__(self, name=None):
        """Creates an instance of the class."""

        self.name = "Unnamed"

        self._model = 0
        self._bonds = None
        self._kdtree = None
        self._atoms = None
        self._coords = None

    # Internal dunder methods
    def __str__(self):
        """String representation of the Structure."""
        return f"<Structure name='{self.name}' natoms={self.natoms}>"

    def __iter__(self):
        """Returns a generator to iterate over the children atoms."""
        yield from self._atoms

    def __len__(self):
        """Point users in the right direction."""
        return (
            "Length of Structure is ambiguous: use natoms or nmodels."
        )

    # Public Methods
    def atom(self, serial):
        """Returns an Atom from the structure.

        Safely translates atom serial to atom array index, handling structures
        with DisorderedAtoms.

        Arguments
        ---------
            serial : int
                serial number of the atom to retrieve.
        """

        try:
            idx = self._idxdict[serial]
        except KeyError:
            max_serial = max(self._idxdict)
            raise StructureError(
                f"Unknown atom: #{serial}. Known indexes: 0 to {max_serial}"
            )

        return self._atoms[idx]

    @property
    def atoms(self):
        """Returns a generator over the Structure's atoms."""
        yield from self._atoms

    @atoms.setter
    def atoms(self, atom_list):
        """Sets atoms as children of this structure.

        Creates an internal mapping between natural indexes (0 .. nth-atom) and
        serial numbers, to allow proper indexing with DisorderedAtoms. In
        addition, this method also clears the coordinate array of the structure,
        the logic being atoms come first, coordinates afterwards.

        Arguments
        ---------
            atom_list : iterable of Atom and/or DisorderedAtom objects
                metadata for each of the atoms in the structure.
        """

        idxdict = {}

        for idx, atom in enumerate(atom_list):
            if isinstance(atom, DisorderedAtom):
                for altloc in atom:
                    idxdict[altloc.serial] = idx
                    altloc.parent = self
            elif isinstance(atom, Atom):
                idxdict[atom.serial] = idx
                atom.parent = self
            else:
                raise StructureError(
                    f'Type of item {idx} is not supported: {type(atom)}'
                )

        self._idxdict = idxdict
        self._atoms = atom_list

    def unpack_atoms(self):
        """Returns a generator of all atoms (including altlocs)."""
        for atom in self:
            if isinstance(atom, DisorderedAtom):
                for altloc in atom:
                    yield altloc
            else:
                yield atom

    @property
    def natoms(self):
        """Returns the number of atoms in the structure"""
        return len(self._atoms)

    @property
    def model(self):
        """Returns the index of the active model."""
        return self._model

    @model.setter
    def model(self, index):
        """Defines which model to pick data from in multi-model structures.

        Arguments
        ---------
            index : int
                (0-based) index of a model to select.
        """

        if not (0 <= index < self.nmodels):
            emsg = f"Model index must be between 0 and {self.nmodels}"
            raise StructureError(emsg)

        self._model = index
        self._kdtree = None  # delete cached kdtree
        logging.debug(f"Set active model to: {index}")

    @property
    def nmodels(self):
        """Returns the number of models in the structure"""
        return self._coords.shape[1]

    @property
    def coords(self):
        """Returns a view of the coordinate array for the selected model."""
        return self._coords[:, self._model]

    @coords.setter
    def coords(self, coord_array):
        """Sets the coordinate array for the entire Structure

        Arguments
        ---------
            coord_array : np.ndarray of np.float
                array of shape NxMx3, where M is the number of models in the
                Structure and N is the number of atoms. The number of atoms
                _must_ match those already pre-loaded in the Structure.
        """

        assert isinstance(coord_array, np.ndarray)
        assert coord_array.ndim == 3, \
            f'Coordinate array has wrong dimensions: {coord_array.ndim} != 3'

        n_total_atoms = len(list(self.unpack_atoms()))  # expand disorders
        assert coord_array.shape[0] == n_total_atoms, \
            f'{coord_array.shape[0]} != {n_total_atoms}'

        assert coord_array.shape[2] == 3, \
            f'{coord_array.shape[2]} != 3'

        self._coords = coord_array

    @property
    def raw_coords(self):
        """Returns the entire coordinate array (all models)."""
        return self._coords

    #
    # Interfaces with external objects
    #
    def _make_kdtree(self):
        """Creates/Returns a KDTree for fast neighbor search.

        KDTree implementation adapted from Biopython by Michiel de Hoon.
        For details, read the source code at interfacea/src/kdtrees.c
        """

        if self._kdtree:
            logging.debug(f"Deleting existing KDTree")
            del self._kdtree

        # Select coordinates of active model/atoms as double-precision floats
        # as required by the C code.
        xyz = np.asarray(self.coords, order='C', dtype=np.float64)
        self._kdtree = kdtrees.KDTree(xyz)

        logging.debug(f"Created KDTree for model #{self.model}")
        del xyz  # no need to keep the copy in memory

    def _neighbors(self, atom, radius):
        """Accesses the KD-tree structure to find neighbors.

        Helper method for self.neighbors.
        """

        dupes = set((atom.serial,))
        for point in self._kdtree.search(atom.coords, radius):

            atom = self.atom(point.index)
            if atom.serial in dupes:
                continue

            # Add serial, not raw index, to avoid multiple altloc matches.
            dupes.add(atom.serial)

            yield (atom, point.radius)

    def neighbors(self, atom, radius):
        """Finds all Atoms within the given radius of itself.

        Arguments
        ---------
            atom : Atom or DisorderedAtom
                central atom to find neighbors of.
            radius : float
                distance cutoff in Angstrom to define neighbors.

        Returns
        -------
            a generator with 2-item tuples, each containing a neighbor Atom and
            its the distance to the query.
        """

        if not isinstance(atom, (Atom, DisorderedAtom)):
            raise StructureError(
                f'First argument "{atom}" is not of type Atom or DisorderedAtom'
            )

        if radius <= 0:
            raise StructureError(
                f"Radius is not a positive, non-zero number: {radius}"
            )

        if self._kdtree is None:
            self._make_kdtree()

        return self._neighbors(atom, radius)

    @property
    def bonds(self):
        """Returns a bond graph describing atom connectivity.

        We use SimpleBondDetector if the structure has no bond graph defined. It
        returns a graph holding atoms as nodes with a Z (atomic number)
        attribute.
        """

        if self._bonds is None:
            ba = SimpleBondDetector(self)
            self.bonds = ba.run()

        return self._bonds

    @bonds.setter
    def bonds(self, bg):
        """Sets a bond graph describing atom connectivity

        Arguments
        ---------
            bg : nx.Graph
                graph describing bonds between atoms in the structure. See
                interfacea.chemistry.bonds for details on the structure of the
                graph.
        """

        assert isinstance(bg, nx.Graph), \
            f"Bond graph object must be a nx.Graph instance: {type(bg)}"

        assert len(bg) == self.natoms, \
            f"Num. of graph nodes != Num. of atoms: {len(bg)} != {self.natoms}"

        self._bonds = bg
        logging.debug(f'Set {bg} as structure bondgraph')
