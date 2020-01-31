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

Structures are meant to be manipulated but _not created_ by users. Refer to the
classes within the io module for that purpose.
"""

import logging
# import re
import weakref

import numpy as np

from ..exceptions import (
    OrphanAtomError
)

logging.getLogger(__name__).addHandler(logging.NullHandler())


class DisorderedAtom(object):
    """Wrapper class for several Atoms sharing the same metadata.

    Allows seamless interaction with a selected instance (altloc) of
    this atom while storing information on the disordered states.
    """

    def __init__(self, parent):

        self.children = {}  # holds Atom objects
        self.selected_child = None

    def __repr__(self):
        """Pretty printing."""
        return f"<DisorderedAtom name={self.name}>"

    def __getattr__(self, attr):
        """Forward calls to selected child."""
        return getattr(self.selected_child, attr)

    def __setattr__(self, attr, value):
        """Forward calls to selected child."""
        return setattr(self.selected_child, attr, value)

    def add(self, atom):
        """Adds an Atom object.

        Args:
            atom (Atom): child object to add.

        Raises:
            NotAnAtomError: if the object being added is not an Atom instance.
            DuplicateAltLocError: if there is already a child with this
                altloc identifier in the DisorderedAtom object.
        """
        if isinstance(atom, Atom):
            if atom.altloc in self.children:
                emsg = f"Altloc '{atom.altloc}' already exists in {self}"
                raise DuplicateAltLocError(emsg) from None
        else:
            emsg = f"{atom} is not an Atom: {type(atom)}"
            raise NotAnAtomError(emsg) from None

        self.children[atom.altloc] = atom

    def select(self, altloc):
        """Selects a child as representative.

        Args:
            altloc (str): altloc identifier of the child Atom.

        Raises
            AtomNotFoundError: if altloc is not in the DisorderedAtom wrapper.
        """
        try:
            self.selected_child = self.children[altloc]
        except KeyError:
            emsg = f"Alternate location '{altloc}' not found in {self}"
            raise AtomNotFoundError(emsg) from None


class Atom(object):
    """Container class to store atomic metadata.

    Args:
        name (str): string to identify the atom (e.g. 'CA').
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

    Attributes:
        coord   (np.array): array of shape (,3) representing cartesian
            coordinates of the atom in Angstrom.
        element (Element): atomic element.
    """

    @classmethod
    def from_atomrecord(cls, record):
        """Creates an Atom class instance from an AtomRecord dataclass.

        Args:
            atomdata (io.AtomRecord): data class containing information to
                create the Atom object
        """

        attrs = record.__dict__
        del attrs['name']
        del attrs['serial']
        return cls(record.name, record.serial, attrs)

    def __init__(self, name, serial, **kwargs):
        """Manually instantiates an Atom class instance."""

        self._parent = None
        self._coord = None

        self.name = name
        self.serial = serial

        self.__dict__.update(kwargs)

        # self._guess_element()

    @property
    def parent(self):
        """Structure to which the Atom belongs to, if any."""
        return self._parent

    @parent.setter
    def parent(self, value):
        if value is not None and isinstance(value, Structure):
            self._parent = weakref.ref(value)  # avoid uncollected garbage

    @property
    def coord(self):
        """Cartesian coordinates of the atom, stored in the parent Structure"""
        if self._parent is not None:
            return self.parent.coord[self.serial]

        emsg = f"Atom does not have a parent Structure with coordinate data"
        raise OrphanAtomError(emsg)

    @coord.setter
    def coord(self, value):
        if value is not None:
            self.coord = value

    # def _guess_element(self):
    #     """Guesses atomic element from atom name."""

    #     alpha = re.sub("[^a-zA-Z]", "", self.name)


class Structure(object):
    """Represents 3D molecules as collections of atoms.

    Args:
        name (str): string that identifies the structure.
        natoms (int): number of atoms in the structure.
        nmodels (int, optional): number of models/frames in the structure.
            Defaults to 1.

    Attributes:
        atoms (list): ordered list of all atoms belonging to the structure.
        coord (np.ndarray): array of shape (nmodels, ndatoms, 3)
    """

    def __init__(self, name, natoms, nmodels=1, **kwargs):
        """Creates an instance of the class."""

        self.name = name
        self.natoms = natoms
        self.nmodels = nmodels

        self.precision = kwargs.get('precision', np.float16)
        self.coords = np.array(
            (self.nmodels, self.natoms, 3),
            dtype=self.precision)

        self.atoms = []
        self._atom_dict = {}  # for __getitem__

    # Internal dunder methods
    def __hash__(self):
        return hash(self.coord.tobytes())  # coordinate hash

    def __repr__(self):
        """String representation of the Structure."""
        return f"<Structure name='{self.name}' natoms={self.natoms}>"

    def __getitem__(self, key):
        """Returns an atom from the structure.

        Args:
            key (str): The key is a unique identifier for an Atom formed by
                four elements: 'chain:residue number:icode:atom name'.
                e.g. 'A:65::CA', 'X:12:B:CB'

        Returns:
            Atom object corresponding to the key.

        Raises:
            AtomNotFoundError: the specified key was not found in the
                Structure object.
        """

        try:
            return self._atom_dict[key]
        except KeyError:
            emsg = f"Atom '{key}' not found in structure"
            raise AtomNotFoundError(emsg) from None

    # 'Private' Methods

    # Public Methods


# Exceptions

class AtomNotFoundError(Exception):
    pass


class DuplicateAltLocError(Exception):
    pass


class NotAnAtomError(Exception):
    pass
