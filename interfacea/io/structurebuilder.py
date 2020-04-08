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
Module containing factory class to build Structures.
"""

import logging

from interfacea.core.atom import Atom, DisorderedAtom
from interfacea.exceptions import StructureBuilderError

logging.getLogger(__name__).addHandler(logging.NullHandler())


class StructureBuilder(object):
    """Factory class for Structure objects.

    This class is used by Readers to build Structure objects.

    Arguments
    ---------
        name : str
            name of the structure being built.

        skip_altloc : bool, optional
            Controls whether alternate locations of the same atom create
            DisorderedAtom objects. If True, builder adds only the first
            instance of each atom. Set to False by default.

    Attributes
    ----------
        atoms : list
            Atom objects to be added to the structure upon calling build().
        coord : list
            3D coordinates of the atoms in the structure.

    """

    def __init__(self, name, **kwargs):

        self.serial = 0  # serial is (re)set automatically from 0

        self.name = name
        self.params = kwargs

        self.atoms = []
        self.coord = []

        self._atoms = {}  # unique id -> serial (to handle DisorderdAtoms)

    def clear(self):
        """Restores the StructureBuilder to its initial state."""

        # Is there a better way?
        self.serial = 0
        self.atoms = []
        self.coord = []
        self._atoms = {}

    def _add_single_atom(self, atom):
        """Adds a single Atom object to the atom list.

        Arguments
        ---------
            atom : Atom
                the atom in question.
        """

        self._atoms[atom.unique_id] = len(self.atoms)
        self.atoms.append(atom)

        self.serial += 1

        logging.debug(f'Added atom to structure: {atom}')

    def _add_disordered_atom(self, idx, atom):
        """Replaces the specified Atom with a DisorderedAtom in the atom list.

        Arguments
        ---------
            idx: int
                the numerical index of the atom in the atom list.

        Returns
        -------
            a DisorderedAtom object with Atom as a child.
        """

        atom_at_idx = self.atoms[idx]
        if isinstance(atom_at_idx, Atom):
            try:
                da = DisorderedAtom()
                da.add(atom_at_idx)
                atom_at_idx = self.atoms[idx] = da

                logging.debug(f'Created DisorderedAtom for atom: {atom_at_idx}')

            except Exception as err:
                raise StructureBuilderError(
                    f'Cannot create a DisorderedAtom from {atom}'
                    f' of type {type(atom)}: {err}'
                )

        atom_at_idx.add(atom)
        logging.debug(f'Added {atom} as a child of {atom_at_idx}')

    def add_atom(self, name, metadata):
        """Adds an atom to the structure.

        Arguments
        ---------
            name : str
                name of the atom, stripped of spaces.
            metadata : dict
                dictionary with properties for the atom, e.g. chain, resid, etc.
                Each property will be translated to an attribute of the Atom
                object. For more information on Atom objects, read the
                documentation on interfacea.core.atom.
        """

        try:
            atom = Atom(name, serial=self.serial, **metadata)
        except Exception as err:
            emsg = f'Failed to create atom #{self.serial}: {name} ({metadata})'
            raise StructureBuilderError(emsg) from err

        idx = self._atoms.get(atom.unique_id)

        if idx is None:  # new atom
            self._add_single_atom(atom)

        elif not self.params.get('skip_altloc'):  # new altloc for existing atom
            self._add_disordered_atom(idx, atom)

        else:  # ignore
            logging.debug(f'Ignored altloc for {atom}')

    def add_coords(self, coords):
        """Adds a coordinate array to the structure.

        The number of coordinates must match the number of atoms already in the
        structure. For multi-model structures, provide a list of lists that must
        be of the same size.

        Arguments
        ---------
            coords : list (of lists) of float
        """

        # is list of lists or list of floats?
        if all(map(lambda x: isinstance(x, list)), coords):  # ensemble
            pass
        if all(map(lambda x: isinstance(x, float)), coords):  # single
            pass

    def build(self):
        """Builds and returns a complete Structure object"""
        raise NotImplementedError('Get to work!')
