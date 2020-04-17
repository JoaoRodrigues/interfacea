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

import numpy as np

from interfacea.io.base import (
    BaseStructureBuilder,
    StructureBuilderError
)

from interfacea.core.atom import Atom, DisorderedAtom
from interfacea.core.structure import Structure

logging.getLogger(__name__).addHandler(logging.NullHandler())


class AllAtomStructureBuilderError(StructureBuilderError):
    pass


class AllAtomStructureBuilder(BaseStructureBuilder):
    """Factory class for all-atom Structures.

    Includes all atoms present in the structure, including disordered atoms. All
    altlocs are converted to DisorderedAtom objects.

    Arguments
    ---------
        parser : BaseParser
            instance of BaseParser subclass.
    """

    def __init__(self, parser):

        self.parser = parser
        self.structure = Structure()

    def _build(self):
        """Validates parsed data and assigns it to the structure object"""

        # Run parser
        self.parser.parse()

        # Assign to structure
        self.structure.atoms = self._add_atoms(self.parser.atoms)
        self.structure.coords = np.asarray(self.parser.coords)

    def _add_atoms(self, parsed_atoms):
        """Adds parsed atoms to the atom list"""

        self._atom_dict = {}
        self._atoms = []

        for atom in parsed_atoms:

            idx = self._atom_dict.get(atom.id)
            if idx is None:  # new atom
                self._add_new_atom(atom)
            else:  # new altloc for an existing atom
                self._add_atom_altloc(idx, atom)

        return self._atoms

    def _add_new_atom(self, atom):
        """Adds a single Atom to the atom list"""

        self._atom_dict[atom.id] = len(self._atoms)
        self._atoms.append(atom)
        logging.debug(f'Atom added to structure: {atom}')

    def _add_atom_altloc(self, idx, atom):
        """Adds a new altloc of an existing Atom to the atom list.

        Makes DisorderedAtom object at self._atoms[idx] if necessary.

        Arguments
        ---------
            idx : int
                internal index of the existing Atom in the atom list.
            atom : Atom
                new atom object to add as a new altloc.
        """

        atom_at_idx = self._atoms[idx]

        # Make DisorderedAtom at idx?
        if isinstance(atom_at_idx, Atom):
            try:
                da = DisorderedAtom()
                da.add(atom_at_idx)
                atom_at_idx = self._atoms[idx] = da

                logging.debug(f'Created DisorderedAtom for atom: {atom_at_idx}')

            except Exception as err:
                raise StructureBuilderError(
                    f'Cannot create a DisorderedAtom from {atom_at_idx}'
                    f' of type {type(atom_at_idx)}: {err}'
                ) from None

        # Add altloc
        try:
            atom_at_idx.add(atom)
            logging.debug(f'Added {atom} as an altloc of {atom_at_idx}')
        except Exception as err:
            raise StructureBuilderError(
                f'Cannot add "{atom}" to disordered atom "{atom_at_idx}": {err}'
            ) from None

    def build(self):
        """Returns a Structure object built from the parser data."""
        self._build()
        return self.structure
