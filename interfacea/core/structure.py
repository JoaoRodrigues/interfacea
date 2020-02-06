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

import collections
import logging
# import re
import weakref

import numpy as np

from interfacea.exceptions import (
    DuplicateAltLocError,
    StructureBuildError,
)

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
    """

    def __init__(self, name, serial, **kwargs):
        """Manually instantiates an Atom class instance."""

        self._parent = None
        self._coords = None

        self.name = name
        self.serial = serial

        self.__dict__.update(kwargs)

    # Dunder methods
    def __str__(self):
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
        return cls(record.name, record.serial, **attrs)

    @property
    def parent(self):
        """Returns the structure the Atom belongs to or None if unbound."""
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
        del self._parent

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
        return f"<Disordered name={self.name} naltlocs={len(self.children)}>"

    def __iter__(self):
        """Returns an iterator over child atoms."""
        for child in self.children.values():
            yield child

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

        self.children[atom.altloc] = atom

        if not (self.selected_child and self.selected_child.occ >= atom.occ):
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


###############################################################################
class Structure(object):
    """Represents 3D molecules as collections of atoms.

    Args:
        name (str): string that identifies the structure.
        coords (np.ndarray): array of shape (num_models, num_atoms, 3)
        atoms (list): Atom objects belonging to the structure, ordered
            by serial number.

    Attributes:
        atoms       (list): ordered list of all atoms in the structure.
        precision   (np.dtype): numerical precision of the atomic coordinates.
        natoms      (int): number of atoms in the structure.
        nmodels     (int): number of models in the structure.
        active_model    (int): 0-based index of the active model.
    """

    def __init__(self, name, coords, atoms):
        """Creates an instance of the class."""

        assert isinstance(coords, np.ndarray) and coords.ndim == 3

        self.name = name
        self._coords = coords
        self.atoms = atoms

        self.precision = coords.dtype
        self.nmodels, self.natoms, _ = coords.shape

        self.model = 0  # default

        self._bind_atoms()
        self._make_atom_dict()

    # Class method to build structure from parser data.
    @classmethod
    def build(cls, name, atomrecords, **kwargs):
        """Creates and returns a Structure object from AtomRecord objects.

        Args:
            name (str): string to use as the resulting Structure name.
            atomrecords (list): list of Atom objects to include in the Structure.

            discard_altloc (bool, optional): ignore atoms with more than one
                instance, i.e. partial occupancies. If True, keeps the altloc
                with the highest occupancy value (and the first in case of a
                tie). If False, builds DisorderedAtom wrappers when needed.
                Default is False.
            precision (np.dtype, optional): numerical precision for storing
                atomic coordinates. Default is np.float16.

        Returns:
            A new Structure object with metadata and coordinate data.

        Raises:
            something
        """

        _args = {
            'precision': np.float16,
            'discard_altloc': False,
        }
        _args.update(kwargs)

        # Pack coordinates
        coords = cls._build_coord_array(
            atomrecords,
            _args.get('precision')
        )

        # Filter altlocs
        record_dict = collections.defaultdict(list)  # uniq -> [loc1, ..]
        for r in atomrecords:
            uniq_id = (r.model, r.chain, r.resid, r.icode, r.name)
            record_dict[uniq_id].append(r)

        if _args.get('discard_altloc'):
            atoms, coords = cls._build_without_altlocs(
                record_dict,
                coords
            )
        else:
            atoms = cls._build_with_altlocs(
                record_dict,
            )

        return cls(name, coords, atoms)

    # Auxiliary methods for build
    @staticmethod
    def _build_coord_array(reclist, dtype):
        """Creates a coordinate array from a list of AtomRecords.

        Args:
            reclist (list): list of AtomRecords with coordinate data.
            dtype (np.dtype): data type for the numpy array.

        Returns:
            A MxAx3 numpy array where M is the number of models and A is the
            number of atoms in the structure.
        """
        coord_list = []  # [ [model1], [model2], ...]
        for rec in reclist:
            if rec.model == len(coord_list):
                coord_list.append([])

            coord_list[-1].append((rec.x, rec.y, rec.z))

        return np.asarray(
            coord_list,
            dtype=dtype
        )

    @staticmethod
    def _build_without_altlocs(recdict, coords):
        """Keeps only highest occupancy instances when altlocs exist.

        Args:
            recdict (dict): dictionary with unique atom identifiers as keys and
                lists of AtomRecord objects as values.
            coords  (np.array): array with coordinate data for all atoms.

        Returns:
            atoms   (list): a list of Atom objects.
            coords  (np.array): a filtered coordinate array without data for
                altlocs.
        """

        atoms = []

        to_remove = set()
        for r in recdict:
            locs = recdict[r]
            if len(locs) > 1:
                locs.sort(key=lambda x: x.occ, reverse=True)
                locs.sort(key=lambda x: x.serial)
                to_remove.update((l.serial for l in locs[1:]))

            atom = Atom.from_atomrecord(locs[0])
            atom.altloc = ''
            atoms.append(atom)

        # Filter coordinates
        coords = np.delete(coords, sorted(to_remove), 1)

        return (atoms, coords)

    @staticmethod
    def _build_with_altlocs(recdict):
        """Builds DisorderedAtom wrappers if necessary.

        Args:
            recdict (dict): dictionary with unique atom identifiers as keys and
                lists of AtomRecord objects as values.

        Returns:
            atoms   (list): a list of Atom/DisorderedAtom objects.
        """

        atoms = []
        for r in recdict:
            locs = recdict[r]
            if len(locs) > 1:
                atom = DisorderedAtom()
                atomlist = (Atom.from_atomrecord(l) for l in locs)
                atom.from_list(atomlist)
            else:
                atom = Atom.from_atomrecord(locs[0])

            atoms.append(atom)

        return atoms

    # Internal dunder methods
    def __str__(self):
        """String representation of the Structure."""
        return f"<Structure name='{self.name}' natoms={self.natoms}>"

    def __iter__(self):
        """Returns a generator to iterate over the children atoms."""
        for atom in self.atoms:
            yield atom

    def __getitem__(self, key):
        """Returns an Atom from the structure.

        Args:
            key (str): The key is a unique identifier for an Atom formed by
                four elements: 'chain:residue number:icode:atom name'.
                e.g. 'A:65::CA', 'X:12:B:CB'

        Returns:
            Atom object corresponding to the input key.

        Raises:
            KeyError: the atom was not found in the structure.
        """

        try:
            return self.atoms[self._atom_dict[key]]
        except KeyError:
            emsg = f"Atom '{key}' not found in structure"
            raise KeyError(emsg) from None

    # 'Private' Methods
    def _bind_atoms(self):
        """Attaches current Atom objects to this Structure object."""

        for atom in self.unpack_atoms():
            atom.parent = self

    def _make_atom_dict(self):
        """Builds a mapping to retrieve individual atoms directly."""

        keys = ('chain', 'resid', 'icode', 'name')  # enough to uniq

        _atom_dict = {}
        for a in self.atoms:
            key = ':'.join(map(str, (getattr(a, k, '') for k in keys)))
            _atom_dict[key] = a.serial

        nuniq = len(_atom_dict)
        if nuniq != self.natoms:
            emsg = (
                "Ambiguous atom identifiers:"
                f"{nuniq} unique vs {self.natoms} total"
            )
            raise StructureBuildError(emsg)

        self._atom_dict = _atom_dict
        logging.debug('Built atom mapping for __getitem__')

    # Public Methods
    def unpack_atoms(self):
        """Returns a generator over all atoms, including all altlocs."""
        for atom in self:
            if isinstance(atom, DisorderedAtom):
                yield from atom
            else:
                yield atom

    @property
    def model(self):
        """Returns the active model."""
        return self._active_model

    @model.setter
    def model(self, index):
        """Defines which model to pick data from in multi-model structures.

        Args:
            index (int): (0-based) index of the model to activate.
        """

        if not (0 <= index < self.nmodels):
            emsg = f"Model index must be between 0 and {self.nmodels}"
            raise ValueError(emsg)

        self._active_model = index

    def _bad_access(self):
        """Stub to tell users how to modify attributes in-place."""
        emsg = (
            "This property cannot be modified in-place.\n"
            "Bind it to a variable first, e.g.:\n"
            "    xyz = s.coords\n"
            "    xyz -= [100, 0, 0]\n"
        )
        raise NotImplementedError(emsg)

    @property
    def coords(self):
        """Returns a view of the coordinate array for the active model."""
        return self._coords[self._active_model]

    @coords.setter
    def coords(self, value):
        self._bad_access()

    @property
    def coords_array(self):
        """Returns a view of the entire coordinate array (all models)."""
        return self._coords

    @coords_array.setter
    def coords_array(self, value):
        self._bad_access()
