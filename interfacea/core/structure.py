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

from interfacea.core.atom import DisorderedAtom
from interfacea.chemistry.bonds import SimpleBondAnalyzer
from interfacea.spatial import kdtrees

logging.getLogger(__name__).addHandler(logging.NullHandler())


class Structure(object):
    """Represents 3D molecules as collections of atoms.

    Arguments
    ---------
        name: str
            string that identifies the structure.
        atoms : list
            Atom objects belonging to the structure, ordered by serial number.
        coords : np.ndarray
            array of shape (num_models, num_atoms, 3)

        make_kdtree : bool, optional
            automatically generate kdtree on instance creation. KDTree objects
            cannot be pickled/copied. Default is True.

    Attributes
    ----------
        atoms : list
            ordered list of all atoms in the structure.
        bonds : networkX.Graph
            bond graph of all atoms in the structure.
        coords : np.ndarray
            coordinate array of the selected model.
        all_coords : np.ndarray
            coordinate array of all models.
        natoms : int
            number of atoms in the structure.
        nmodels : int
            number of models in the structure.
    """

    def __init__(self, name, atoms, coords, **kwargs):
        """Creates an instance of the class."""

        _defaults = {
            'make_kdtree': True,
        }

        self.__dict__.update(_defaults)
        self.__dict__.update(kwargs)

        # Sanity check if someone decides to make their own Structure
        assert isinstance(coords, np.ndarray) and coords.ndim == 3

        self.name = name
        self._add_atoms(atoms)

        self._coords = coords
        self._model = 0

        self._bonds = None
        self._kdtree = None

        self._make_kdtree()

    # Class method to build structure from parser data.
    # @classmethod
    # def from_atomrecords(cls, name, atomrecords, **kwargs):
    #     """Creates and returns a Structure object from AtomRecord objects.

    #     Arguments
    #     ---------
    #         name : str
    #             name of the resulting structure.
    #         atomrecords : list
    #             AtomRecords to build Atoms from.

    #         discard_altloc : bool, optional
    #             ignore atoms with more than one alternate locations. Keeps only
    #             the first altloc. If False (default), builds DisorderedAtom
    #             wrappers to store multiple locations.

    #     Returns
    #     -------
    #         Structure instance with metadata and coordinate data.
    #     """

    #     _defaults = {
    #         'discard_altloc': False,
    #     }
    #     _defaults.update(kwargs)

    #     atoms, coords = [], []
    #     record_dict = {}  # uniq -> serial
    #     for r in atomrecords:
    #         uniq_id = (r.model, r.chain, r.resid, r.icode, r.name)
    #         idx = record_dict.get(uniq_id)

    #         if idx is None:  # new atom
    #             atom = Atom.from_atomrecord(r)
    #             record_dict[uniq_id] = len(atoms)
    #             atoms.append(atom)

    #         elif not _defaults['discard_altloc']:  # new altloc for existing atom
    #             existing = atoms[idx]
    #             if isinstance(existing, Atom):
    #                 logging.debug(f"New disordered atom at #{r.serial}")
    #                 disatom = DisorderedAtom()
    #                 disatom.add(existing)
    #                 atoms[idx] = disatom
    #             new_loc = Atom.from_atomrecord(r)
    #             atoms[idx].add(new_loc)

    #         else:  # ignore
    #             logging.debug(f"Ignoring duplicate atom: #{r.serial}")
    #             continue

    #         if r.model == len(coords):  # make new models
    #             coords.append([])

    #         coords[-1].append((r.x, r.y, r.z))

    #     # Pack coordinates into numpy array
    #     coords = np.asarray(
    #         coords,
    #         dtype=np.float64
    #     )

    #     assert coords.ndim == 3, \
    #         f'Wrong shape for coordinate array: {coords.shape}'

    #     logging.debug(f"Built new structure with {len(atoms)} atoms")
    #     return cls(name, atoms, coords)

    # Internal dunder methods
    def __str__(self):
        """String representation of the Structure."""

        return f"<Structure name='{self.name}' natoms={self.natoms}>"

    def __iter__(self):
        """Returns a generator to iterate over the children atoms."""

        yield from self.atoms

    def __len__(self):
        """Point users in the right direction."""
        return (
            "Length of Structure is ambiguous: use num_atoms or num_models."
        )

    # 'Private' Methods
    def _add_atoms(self, atoms):
        """Sets atoms as childs of this structure.

        Additionally, creates idxdict mapping to avoid IndexErrors when we
        have DisorderedAtoms and try mapping by serial number, ie when
        max(atom.serial) > len(atoms).
        """

        self.atoms = atoms
        self._idxdict = idxdict = {}
        for idx, atom in enumerate(atoms):
            if isinstance(atom, DisorderedAtom):
                for altloc in atom:
                    idxdict[altloc.serial] = idx
                    altloc.parent = self
            else:
                idxdict[atom.serial] = idx
                atom.parent = self

    def _make_kdtree(self):
        """Creates/Returns a KDTree for fast neighbor search.

        KDTree implementation adapted from Biopython by Michiel de Hoon.
        For details, read the source code at interfacea/src/kdtrees.c
        """

        if self._kdtree:
            logging.debug(f"Deleting existing KDTree")
            del self._kdtree

        xyz = self.coords  # only for active model
        xyz = np.asarray(xyz, np.float64)  # kdtrees requires double precision
        self._kdtree = kdtrees.KDTree(xyz)
        logging.debug(f"Created KDTree for model #{self.model}")

    # Public Methods
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
            an iterator with 2-item tuples, each containing a neighbor Atom and
            its the distance to the query.

        Raises
        ------
            ValueError
                when 'radius' is not a positive, non-zero, number.
        """

        if radius <= 0:
            raise ValueError(
                f"Radius is not a positive, non-zero number: {radius}"
            )

        if self._kdtree is None:
            self._make_kdtree()

        kdt = self._kdtree
        coords = np.asarray(atom.coords, dtype=np.float64)  # as 64-bit vector

        dupes = set((atom.serial,))
        for point in kdt.search(coords, radius):

            atom = self.atom(point.index)
            if atom.serial in dupes:
                continue

            # Add serial, not raw index, to avoid multiple altloc matches.
            dupes.add(atom.serial)

            yield (atom, point.radius)

    def unpack_atoms(self):
        """Returns a generator over all atoms, including all altlocs."""

        for atom in self:
            if isinstance(atom, DisorderedAtom):
                yield from atom
            else:
                yield atom

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
            raise KeyError(
                f"Unknown atom: #{serial}. Known indexes: 0 to {max_serial}"
            )

        return self.atoms[idx]

    @property
    def num_atoms(self):
        """Returns the number of atoms in the structure"""
        return len(self.atoms)

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
            raise ValueError(emsg)

        self._model = index
        self._kdtree = None  # delete cached kdtree
        logging.debug(f"Set active model to: {index}")

    @property
    def num_models(self):
        """Returns the number of models in the structure"""
        return self._coords.shape[0]

    @property
    def models(self):
        """Iterates over all models of the structure"""

        start = self.model
        for midx in self.num_models:
            self.model(midx)
            yield self

        self.model(start)  # reset

    @property
    def bonds(self):
        """Returns a bond graph describing atom connectivity.

        We use SimpleBondAnalyzer if the structure has no bond graph defined. It
        returns a graph holding atoms as nodes with a Z (atomic number)
        attribute.
        """

        if self._bonds is None:
            ba = SimpleBondAnalyzer(self)
            self._bonds = ba.run()

        return self._bonds

    @bonds.setter
    def bonds(self, bondgraph):
        """Sets a bond graph describing atom connectivity

        Arguments
        ---------
            bondgraph : nx.Graph
                graph describing bonds between atoms in the structure. See
                interfacea.chemistry.bonds for details on the structure of the
                graph.
        """

        # Check bond graph (roughly) matches structure
        if isinstance(bondgraph, nx.Graph):
            assert len(bondgraph) == self.num_atoms, \
                f"Number of nodes in bond graph differs from number of atoms."

            self._bonds = bondgraph

        raise TypeError(
            f"Bond graph object must be a networkx Graph instance."
        )

    @property
    def coords(self):
        """Returns a view of the coordinate array for the selected model."""
        return self._coords[self._model]

    @property
    def all_coords(self):
        """Returns a view of the entire coordinate array (all models)."""
        return self._coords
