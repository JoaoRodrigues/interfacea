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
            atom = Atom(name, **metadata)
        except Exception as err:
            emsg = f'Failed to create atom #{self.serial}: {name} ({metadata})'
            raise StructureBuilderError(emsg) from err

        atom.serial = self.serial

        # unique id
        atom_uid = tuple(
            getattr(atom, attr, None)
            for attr in ('name', 'chain', 'resid', 'icode')
        )

        idx = self._atoms.get(atom_uid)

        if idx is None:  # new atom
            idx = len(self.atoms)
            self._atoms[atom_uid] = idx
            self.atoms.append(atom)
            logging.debug(f'Added new atom: {atom}')

        elif not self.params.get('skip_altloc'):  # new altloc for existing atom
            previous = self.atoms[idx]
            if isinstance(previous, Atom):  # make new DisorderedAtom
                try:
                    da = DisorderedAtom()
                    da.add(previous)
                except Exception as err:
                    emsg = f'Failed to create disordered atom for {previous}'
                    raise StructureBuilderError(emsg) from err

                self.atoms[idx] = da
                logging.debug(f'Created DisorderedAtom for atom: {previous}')

            da = self.atoms[idx]
            try:
                da.add(atom)
            except Exception as err:
                emsg = f'Failed to add {atom} to {da}'
                raise StructureBuilderError(emsg) from err

            logging.debug(f'Added {atom} as altloc to {da}')

        else:  # ignore
            logging.debug(f'Ignored altloc for {atom}')
            return

        self.serial += 1

    def build(self):
        """Builds and returns a complete Structure object"""
        raise NotImplementedError('Get to work!')
