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
Module containing generic classes to build molecular structures.

The StructureBuilder class should be used by file parsers to build
Structure objects. Centralizes all checking/validation in one place.
"""

import logging

from ..core.Structure import (
    Atom,
    DisorderedAtom,
    Structure,
)

logging.getLogger(__name__).addHandler(logging.NullHandler())


class StructureBuilder(object):
    """Builds a Structure object from data from a parser.
    """

    @classmethod
    def build(cls, name, atomrecords, discard_altloc=False, **kwargs):
        """Loops over the atom data and returns a Structure object.

        Args:
            name (str): string to use as the resulting Structure name.
            atom_data (list): list of Atom objects to include in the Structure.
            discard_altloc (bool, optional): keep only the instance of each atom with
                the highest occupancy value if True. Default is False.

        Returns:
            A new Structure object containing all the metadata but no coordinate data.

        Raises:
            something
        """

        if discard_altloc:
            atomdict = {}
            for atom in atomrecords:
                uniq_id = f"{atom.model}:{atom.chain}:{atom.resid}:{atom.name}"
                existing = atomdict.get(uniq_id)
                if existing and atom.occ <= existing.occ:
                    continue

                atom.altloc = ''
                atomdict[uniq_id] = atom

            atomrecords = sorted(atomdict.values(), key=lambda a: a.serial)
            logging.debug(f"Discarded {len(atomrecords) - len(atomdict)} altlocs")

            nmodels = len({a.model for a in atomrecords})
            atoms = [Atom.from_atomrecord(r) for r in atomrecords]

            s = Structure(name=name, natoms=len(atoms), nmodels=nmodels, kwargs)
            s.coord = np.asarray(
                ((a.x, a.y, a.z) for a in atomrecords),
                dtype=kwargs.get('precision', np.float16)
            )

            return s

        else:  # build DisorderedAtoms if necessary
            raise NotImplementedError
