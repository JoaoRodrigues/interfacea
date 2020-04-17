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

from interfacea.exceptions import BaseInterfaceaException

logging.getLogger(__name__).addHandler(logging.NullHandler())


class DisorderedAtomError(BaseInterfaceaException):
    pass


class Atom(object):
    """Container class to store atomic metadata.

    Arguments
    ---------
        name : str
            string to identify the atom.
        serial : int, optional
            numerical index of the atom in the Structure. Defaults to 0.

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
            array of shape (,3) wih the 3D cartesian coordinates in Angstrom.
            Only defined when bound to a parent Structure.
    """

    def __init__(self, name, serial=0, **kwargs):
        """Manually instantiates an Atom class instance."""

        self._id = None
        self._fid = None
        self._parent = None

        self.name = name
        self.serial = serial

        self.__dict__.update(kwargs)

    # Dunder methods
    def __str__(self):
        """Pretty string representation of the Atom object."""
        return f"<Atom name={self.name} serial={self.serial}>"

    # Public Methods/Attributes
    @property
    def id(self):
        """Returns the id tuple of the atom.

        The id is composed of the atom name, residue number and icode,
        and chain id. Alternate locations of the same atom share the same
        id.

        """
        if self._id is None:
            aid = tuple(
                getattr(self, attr, None)
                for attr in ('name', 'chain', 'resid', 'icode')
            )
            self._id = aid

        return self._id

    @property
    def full_id(self):
        """Returns the full id tuple of the atom.

        The full id is composed of the atom name and altloc, residue number and
        icode, and chain id. Unlike atom.id, this is unique to even altlocs.

        """
        if self._fid is None:
            afid = tuple(
                getattr(self, attr, None)
                for attr in ('name', 'chain', 'resid', 'icode', 'altloc')
            )
            self._fid = afid

        return self._fid

    @property
    def parent(self):
        """Returns the structure the Atom belongs to or None if unbound."""
        try:
            return self._parent()
        except TypeError:  # not callable, not weakref?
            return self._parent

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
    """Class representing a disordered atom, with multiple alternate locations.

    Wraps several Atom objects, allowing transparent interaction with a selected
    child.
    """

    def __init__(self):

        self.__dict__['children'] = {}  # holds Atom objects
        self.__dict__['representative'] = None

    def __str__(self):
        """Pretty printing."""
        if 'name' in self.__dict__:
            return (
                f"<DisorderedAtom name={self.name} "
                f"serial={self.serial} nlocs={self.nlocs}>"
            )

        return f"<Empty DisorderedAtom>"

    def __iter__(self):
        """Returns an iterator over child atoms."""
        yield from self.children.values()

    def __getattr__(self, attr):
        """Forward all unknown calls to selected child."""

        try:
            return getattr(self.representative, attr)
        except AttributeError as err:
            if self.representative is None:
                emsg = "DisorderedAtom has no children."
                raise DisorderedAtomError(emsg)
            raise err from None  # re-raise

    def __setattr__(self, attr, value):
        """Forward all unknown calls to selected child."""
        if attr in self.__dict__:
            self.__dict__[attr] = value
        else:
            setattr(self.representative, attr, value)

    def add(self, atom):
        """Adds an Atom object.

        Arguments
        ---------
            atom : Atom
                child Atom to add to container.
        """

        try:
            altloc = atom.altloc
        except AttributeError:
            altloc = str(self.nlocs)
            atom.altloc = altloc
            logging.warning(
                f'Auto-setting missing altloc for {atom}: {altloc}'
            )

        if altloc in self.children:
            emsg = f"Altloc '{altloc}' already exists in {self}"
            raise DisorderedAtomError(emsg) from None

        self.children[altloc] = atom

        if not self.representative:
            self.representative = atom
            return

        # Replace if atom.occ is larger
        try:
            if self.representative.occ < atom.occ:
                self.representative = atom
        except AttributeError:  # no occ in child or atom
            self.representative = atom
            # see https://docs.python.org/3/howto/logging.html
            logging.warning(f'Could not sort altloc by occupancy for {self}')

    def delete(self, altloc):
        """Removes a child Atom from the DisorderedAtom.

        Arguments
        ---------
            altloc : str
                altloc identifier for the Atom to remove.
        """

        try:
            if self.representative.altloc == altloc:
                self.representative = None  # remove reference
            del self.children[altloc]
        except KeyError:
            raise DisorderedAtomError(
                f'Altloc "{altloc}" not found in DisorderedAtom'
            )
        else:
            if self.children and self.representative is None:
                self.select()

    def from_list(self, atomlist):
        """Adds Atom objects from a list.

        Argument
        --------
            atomlist : list
                obviously, a list of Atoms.
        """
        for atom in atomlist:
            self.add(atom)

    def select(self, altloc=None):
        """Selects a child as representative.

        Argument
        --------
            altloc : str, optional
                altloc identifier of the child Atom. If None, picks the child
                with the highest occupancy. In case of multiple children
                with the same occupancy, picks the first based on alphabetical
                altloc order.
        """

        if not self.children:
            raise DisorderedAtomError('DisorderedAtom has no children')

        if altloc is None:
            self.representative = sorted(
                list(self),
                key=lambda a: (a.occupancy, a.altloc)
            )[0]
        else:
            try:
                self.representative = self.children[altloc]
            except KeyError:
                emsg = f"Alternate location '{altloc}' not found in {self}"
                raise DisorderedAtomError(emsg) from None

        logging.debug(f'Set {self.representative.altloc} as representative')

    @property
    def nlocs(self):
        """Returns the number of children in the DisorderedAtom"""
        return len(self.children)
