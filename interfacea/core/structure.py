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
import weakref

import networkx as nx
import numpy as np

from interfacea.exceptions import (
    DuplicateAltLocError,
)

from interfacea.chemistry.bonds import SimpleBondAnalyzer
from interfacea.spatial import kdtrees

logging.getLogger(__name__).addHandler(logging.NullHandler())


###############################################################################
class Atom(object):
    """Container class to store atomic metadata.

    Args:
        name (str): string to identify the atom.
        serial  (int): numerical index of the atom.

        hetatm  (bool, optional): flag to identify HETATMs.
        altloc  (str, optional): identifier for alternate location.
        resname (str, optional): name of the parent residue.
        resid   (int, optional): residue sequence number.
        icode   (str, optional): residue insertion code.
        chain   (str, optional): name of the parent chain.
        b       (float, optional): temperature factor.
        occ     (float, optional): fractional occupancy.
        segid   (str, optional): segment identifier.
        element (Element, optional): atomic element.

    Attributes:
        coords  (np.array): array of shape (,3) representing cartesian
            coordinates of the atom in Angstrom. Only defined when bound
            to a parent Structure, otherwise raises an error on access.
        is_disordered (bool): True if the atom has multiple locations, i.e.,
            is part of a DisorderedAtom object.
    """

    def __init__(self, name, serial, **kwargs):
        """Manually instantiates an Atom class instance."""

        self._parent = None
        self._coords = None

        self.name = name
        self.serial = serial
        self.is_disordered = False

        self.__dict__.update(kwargs)

    # Dunder methods
    def __str__(self):
        """Pretty string representation of the Atom object."""
        return f"<Atom name={self.name} serial={self.serial}>"

    def __repr__(self):
        """Pretty string representation of the Atom object."""
        return f"<Atom name={self.name} serial={self.serial}>"

    # Public Methods/Attributes
    @classmethod
    def from_atomrecord(cls, record):
        """Creates an Atom class instance from an AtomRecord dataclass.

        Args:
            atomdata (io.AtomRecord): data class containing information to
                create the Atom object
        """

        attrs = record.__dict__.copy()
        del attrs['name']
        del attrs['serial']

        # Ignore x, y, z fields if present
        for f in ('x', 'y', 'z'):
            try:
                del attrs[f]
            except Exception:
                pass

        return cls(record.name, record.serial, **attrs)

    # Public Methods
    # Properties
    @property
    def parent(self):
        """Returns the structure the Atom belongs to or None if unbound."""
        if self._parent is None:
            return None  # unbound
        return self._parent()

    @parent.setter
    def parent(self, value):
        if isinstance(value, Structure):
            self._parent = weakref.ref(value)  # avoid uncollected garbage
        else:
            emsg = f"Parent object must be a Structure type."
            raise TypeError(emsg)

    @parent.deleter
    def parent(self):
        self._parent = None

    @property
    def coords(self):
        """Cartesian coordinates of the atom.

        Raises:
            AttributeError: atom is not bound to a Structure object.
        """
        try:
            return self.parent.coords[self.serial]
        except AttributeError:
            emsg = f"Atom is not bound to a parent Structure."
            raise AttributeError(emsg) from None


###############################################################################
class DisorderedAtom(object):
    """Wrapper class for several Atoms sharing the same metadata.

    Allows seamless interaction with a selected instance (altloc) of
    this atom while storing information on the disordered states.
    """

    def __init__(self):

        self.__dict__['children'] = {}  # holds Atom objects
        self.__dict__['selected_child'] = None

    def __str__(self):
        """Pretty printing."""
        return (
            f"<DisorderedAtom name={self.name} "
            f"serial={self.serial} nlocs={self.nlocs}>"
        )

    def __repr__(self):
        """Pretty printing."""
        return (
            f"<DisorderedAtom name={self.name} "
            f"serial={self.serial} nlocs={self.nlocs}>"
        )

    def __iter__(self):
        """Returns an iterator over child atoms."""
        yield from self.children.values()

    def __getattr__(self, attr):
        """Forward all unknown calls to selected child."""
        if attr in self.__dict__:
            return self.__dict__[attr]

        try:
            return getattr(self.selected_child, attr)
        except AttributeError as err:
            if self.selected_child is None:
                emsg = "DisorderedAtom has no children."
                raise AttributeError(emsg)
            raise err from None  # re-raise

    def __setattr__(self, attr, value):
        """Forward all unknown calls to selected child."""
        if attr in self.__dict__:
            self.__dict__[attr] = value
            return

        setattr(self.selected_child, attr, value)

    def add(self, atom):
        """Adds an Atom object.

        Args:
            atom (Atom): child object to add.

        Raises:
            TypeError: if the object being added is not an Atom instance.
            DuplicateAltLocError: if there is already a child with this
                altloc identifier in the DisorderedAtom object.
        """

        if atom.altloc in self.children:
            emsg = f"Altloc '{atom.altloc}' already exists in {self}"
            raise DuplicateAltLocError(emsg) from None

        atom.is_disordered = True  # flag
        self.children[atom.altloc] = atom

        if not self.selected_child:
            self.selected_child = atom

        # Replace if atom.occ is larger
        if self.selected_child.occ < atom.occ:
            self.selected_child = atom

    def from_list(self, atomlist):
        """Adds Atom objects from a list.

        Args:
            atomlist (list): obviously, a list of Atoms.
        """
        for atom in atomlist:
            self.add(atom)

    def select(self, altloc):
        """Selects a child as representative.

        Args:
            altloc (str): altloc identifier of the child Atom.

        Raises
            KeyError: if altloc is not in the DisorderedAtom wrapper.
        """
        try:
            self.selected_child = self.children[altloc]
        except KeyError:
            emsg = f"Alternate location '{altloc}' not found in {self}"
            raise KeyError(emsg) from None

    @property
    def nlocs(self):
        """Returns the number of children in the DisorderedAtom"""
        return len(self.children)


###############################################################################
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
        self.atoms = atoms

        self._coords = coords
        self._model = 0

        self._kdtree = None
        self._make_kdtree()

        # Bind atoms to structure and map indices
        self._bind_atoms()
        self._map_indices()

    # Class method to build structure from parser data.
    @classmethod
    def from_atomrecords(cls, name, atomrecords, **kwargs):
        """Creates and returns a Structure object from AtomRecord objects.

        Arguments
        ---------
            name : str
                name of the resulting structure.
            atomrecords : list
                AtomRecords to build Atoms from.

            discard_altloc : bool, optional
                ignore atoms with more than one alternate locations. Keeps only
                the first altloc. If False (default), builds DisorderedAtom
                wrappers to store multiple locations.

        Returns
        -------
            Structure instance with metadata and coordinate data.
        """

        _defaults = {
            'discard_altloc': False,
        }
        _defaults.update(kwargs)

        atoms, coords = [], []
        record_dict = {}  # uniq -> serial
        for r in atomrecords:
            uniq_id = (r.model, r.chain, r.resid, r.icode, r.name)
            idx = record_dict.get(uniq_id)

            if idx is None:  # new atom
                atom = Atom.from_atomrecord(r)
                record_dict[uniq_id] = len(atoms)
                atoms.append(atom)

            elif not _defaults['discard_altloc']:  # new altloc for existing atom
                existing = atoms[idx]
                if isinstance(existing, Atom):
                    logging.debug(f"New disordered atom at #{r.serial}")
                    disatom = DisorderedAtom()
                    disatom.add(existing)
                    atoms[idx] = disatom
                new_loc = Atom.from_atomrecord(r)
                atoms[idx].add(new_loc)

            else:  # ignore
                logging.debug(f"Ignoring duplicate atom: #{r.serial}")
                continue

            if r.model == len(coords):  # make new models
                coords.append([])

            coords[-1].append((r.x, r.y, r.z))

        # Pack coordinates into numpy array
        coords = np.asarray(
            coords,
            dtype=np.float64
        )

        assert coords.ndim == 3, \
            f'Wrong shape for coordinate array: {coords.shape}'

        logging.debug(f"Built new structure with {len(atoms)} atoms")
        return cls(name, atoms, coords)

    # Internal dunder methods
    def __str__(self):
        """String representation of the Structure."""

        return f"<Structure name='{self.name}' natoms={self.natoms}>"

    def __iter__(self):
        """Returns a generator to iterate over the children atoms."""

        for atom in self.atoms:
            yield atom

    # 'Private' Methods
    def _bind_atoms(self):
        """Attaches current Atom objects to this Structure object."""

        for atom in self.unpack_atoms():
            atom.parent = self

    def _map_indices(self):
        """Builds a mapping of raw indices to Atom/DisorderedAtom serials."""

        idxdict = {}
        for raw_idx, atom in enumerate(self):
            if isinstance(atom, DisorderedAtom):
                for loc in atom:
                    idxdict[loc.serial] = raw_idx
            else:
                idxdict[atom.serial] = raw_idx
        self._idxdict = idxdict

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

        Args:
            atom (Atom or DisorderedAtom): central atom to find neighbors of.
            radius (float): distance cutoff in Angstrom to define neighbors.

        Returns:
            a generator with 2-item tuples containing the neighboring
            Atom object and its distance to the query.

        Raises:
            TypeError: if 'atom' is not of type Atom or DisorderedAtom.
            ValueError: if 'radius' is not a positive, non-zero, number.
        """

        if not isinstance(atom, (Atom, DisorderedAtom)):
            raise TypeError(
                f"Atom must be of type 'Atom' or 'DisorderedAtom': {type(atom)}"
            )

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

            if point.index not in dupes:
                dupes.add(point.index)
                yield (
                    self.atom(point.index),
                    point.radius
                )

    def unpack_atoms(self):
        """Returns a generator over all atoms, including all altlocs."""
        for atom in self:
            if isinstance(atom, DisorderedAtom):
                yield from atom
            else:
                yield atom

    def atom(self, index):
        """Returns an Atom from the structure.

        Better than accessing self.atoms directly, as it accounts for
        DisorderedAtoms.

        Args:
            index (int): positive integer
        """

        try:
            idx = self._idxdict[index]
        except KeyError:
            raise KeyError(
                f"Unknown atom: #{index}. Known indexes: 0 to {self.num_atoms}"
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
    def coords(self):
        """Returns a view of the coordinate array for the selected model."""
        return self._coords[self._model]

    @property
    def full_coords(self):
        """Returns a view of the entire coordinate array (all models)."""
        return self._coords
