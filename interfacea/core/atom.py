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
Module containing classes to represent atomic metadata in topologies.

Atom classes stored metadata only (names, etc), no coordinates.
"""

import logging

from interfacea.exceptions import InterfaceaError

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class DisorderedAtomError(InterfaceaError):
    """Exception subclass for DisorderedAtom errors."""

    pass


class Atom:
    """Container class to store atomic metadata.

    Arguments
    ---------
        name : str
            string to identify the atom.
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
        occupancy : float, optional
            fractional occupancy.
        segid : str, optional
            segment identifier.
        element : Element, optional
            atomic element.

    Attributes
    ----------
        index : int, optional
            numerical index of the atom if part of a Topology. Default is None.
        id : tuple
            identifier for the atom that includes chain and residue
            information, if existing.
        full_id : tuple
            identifier that includes the altloc field if existing. Should
            uniquely identify the atom.
    """

    def __init__(self, name, **kwargs):
        """Create a new Atom instance."""

        self.name = name
        self.__dict__.update(kwargs)

        self.index = None

        self._id = None
        self._full_id = None

    # Dunder methods
    def __str__(self):
        """Pretty string representation of the Atom object."""

        return f"<Atom {self.name} index={self.index}>"

    def __repr__(self):
        """Pretty printing as well."""

        return self.__str__()

    def __eq__(self, other):
        """Equality between Atom objects.

        DisorderedAtoms will be compared based on their selected child.
        """

        if not isinstance(other, (Atom, DisorderedAtom)):
            return False

        return self.full_id == other.full_id

    # Public Methods/Attributes
    @property
    def id(self):
        """Return the id tuple of the atom."""

        if self._id is None:
            aid = tuple(
                getattr(self, attr, None)
                for attr in ("name", "chain", "resid", "icode")
            )
            self._id = aid

        return self._id

    @property
    def full_id(self):
        """Return the full id tuple of the atom."""
        if self._full_id is None:
            altloc = getattr(self, "altloc", None)
            self._full_id = self.id + (altloc,)

        return self._full_id


class DisorderedAtom(object):
    """Class to represent multiple locations of the same atom.

    Acts as a wrapper for multiple Atom objects, forwarding calls to a selected
    child and behaving as a regular Atom.

    Attributes
    ----------
        children : list
            list of children Atom objects.
        selected : Atom
            primary Atom child to which all methods call will be forwarded.
    """

    # The explicit use of self.__dict__ is to avoid issues with __getattr__ and
    # the forwarding to the selected child.

    def __init__(self):
        """Create a new instance of a DisorderedAtom."""

        self.__dict__["children"] = {}  # maps altlocs to Atom objects
        self.__dict__["selected"] = None

    def __str__(self):
        """Pretty printing."""

        n = len(self)
        if n:
            return f"<DisorderedAtom name={self.name} index={self.index} [{n}]>"
        return "<DisorderedAtom [Empty]>"

    def __repr__(self):
        """Pretty printing as well."""
        return self.__str__()

    def __iter__(self):
        """Return an iterator over child atoms."""

        yield from self.children.values()

    def __len__(self):
        """Return the number of children."""
        return len(self.children)

    def __getattr__(self, attr):
        """Forward all unknown calls to selected child."""

        try:
            return getattr(self.selected, attr)
        except AttributeError as err:
            if self.selected is None:
                raise DisorderedAtomError("DisorderedAtom has no children.")
            raise err from None  # re-raise

    def __setattr__(self, attr, value):
        """Forward all unknown calls to selected child."""
        if attr in self.__dict__:
            self.__dict__[attr] = value
        else:
            setattr(self.selected, attr, value)

    def __getitem__(self, key):
        """Retrieve a child by its altloc identifier."""
        return self.children[key]

    def __eq__(self, other):
        """Return True if two DisorderedAtoms are equivalent.

        When comparing an Atom to a DisorderedAtom, compares the Atom with the
        selected child.
        """

        if isinstance(other, Atom):  # compare to top child
            return self.selected == other
        elif isinstance(other, DisorderedAtom):
            if len(self) != len(other):
                return False

            for a1, a2 in zip(self.children, other.children):
                if a1 != a2:
                    return False

            return True
        else:
            return False

    def add(self, atom):
        """Add a child Atom object."""

        try:
            altloc = atom.altloc
        except AttributeError:
            altloc = str(len(self))
            atom.altloc = altloc
            logger.warning(f"Auto-setting missing altloc for {atom}: {altloc}")

        if altloc in self.children:
            emsg = f"Altloc '{altloc}' already exists in {self}"
            raise DisorderedAtomError(emsg) from None

        self.children[altloc] = atom

        if not self.selected:
            self.selected = atom
            return

        # Replace if atom.occupancy is larger
        try:
            if self.selected.occupancy < atom.occupancy:
                self.selected = atom
        except AttributeError:  # no occupancy in child or atom
            self.selected = atom
            # see https://docs.python.org/3/howto/logging.html
            logger.warning(f"Could not sort altloc by occupancy for {self}")

    def delete(self, altloc):
        """Remove a child Atom from the DisorderedAtom.

        Arguments
        ---------
            altloc : str
                altloc identifier for the Atom to remove.
        """

        try:
            if self.selected.altloc == altloc:
                self.selected = None  # remove reference
            del self.children[altloc]
        except KeyError:
            raise DisorderedAtomError(f'Altloc "{altloc}" not found in DisorderedAtom')
        else:
            if self.children and self.selected is None:
                self.select()

    def from_list(self, atomlist):
        """Add Atom objects from a list.

        Arguments
        ---------
            atomlist : list
                obviously, a list of Atoms.
        """
        for atom in atomlist:
            self.add(atom)

    def select(self, altloc=None):
        """Define a child as the selected child of the DisorderedAtom.

        Arguments
        ---------
            altloc : str, optional
                altloc identifier of the child Atom. If None, picks the child
                with the highest occupancy. In case of multiple children
                with the same occupancy, picks the first based on alphabetical
                altloc order.
        """

        if not self.children:
            raise DisorderedAtomError("DisorderedAtom has no children")

        if altloc is None:
            self.selected = sorted(self, key=lambda a: (a.occupancy, a.altloc))[-1]
        else:
            try:
                self.selected = self.children[altloc]
            except KeyError:
                emsg = f"Alternate location '{altloc}' not found in {self}"
                raise DisorderedAtomError(emsg) from None

        logger.debug(f"Set {self.selected.altloc} as selected")
