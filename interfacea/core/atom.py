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
Module containing classes to represent atomic metadata.

Atom and related classes stored metadata only (names, serial numbers, etc),
since coordinate information is kept at the Structure level.
"""

import logging
import weakref

from interfacea.exceptions import (
    DuplicateAltLocError,
)

logging.getLogger(__name__).addHandler(logging.NullHandler())


class Atom(object):
    """Container class to store atomic metadata.

    Arguments
    ---------
        name : str
            string to identify the atom.
        serial : int
            numerical index of the atom.

        hetatm : bool, optional
            flag to identify HETATMs.
        altloc : str, optional
            identifier for alternate location.
        resname : str, optional
            name of the parent residue.
        resid : int, optional
            residue sequence number.
        icode : str, optional
            residue insertion code.
        chain : str, optional
            name of the parent chain.
        b : float, optional
            temperature factor.
        occ : float, optional
            fractional occupancy.
        segid : str, optional
            segment identifier.
        element : Element, optional
            atomic element.

    Attributes
    ----------
        coords : np.array
            array of shape (,3) representing cartesian coordinates of the atom
            in Angstrom. Only defined when bound to a parent Structure,
            otherwise raises an error on access.
        is_disordered : bool
            True if the atom has multiple locations, i.e., is part of a
            DisorderedAtom object.
    """

    def __init__(self, name, serial, **kwargs):
        """Manually instantiates an Atom class instance."""

        self._parent = None

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

        Argument
        --------
            atomdata : io.AtomRecord
                data class containing information to create the Atom object.
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
        self._parent = weakref.ref(value)  # avoid uncollected garbage

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
        except AttributeError as err:
            raise AttributeError(f"Is this atom bound to a Structure?") from err


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

        Arguments
        ---------
            atom : Atom)
                child Atom to add to container.

        Raises
        ------
            TypeError
                when the object being added is not an Atom instance.
            DuplicateAltLocError
                when there is already a child with this altloc identifier in the
                DisorderedAtom object.
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

        Argument
        --------
            atomlist : list
                obviously, a list of Atoms.
        """
        for atom in atomlist:
            self.add(atom)

    def select(self, altloc):
        """Selects a child as representative.

        Argument
        --------
            altloc : str
                altloc identifier of the child Atom.

        Raises
        ------
            KeyError
                when altloc is not in the DisorderedAtom wrapper.
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
