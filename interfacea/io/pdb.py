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
Module containing classes to build molecular structures from PDB files.
"""

import collections
import logging
import warnings

from interfacea.exceptions import (
    PDBReaderError,
    PDBReaderWarning
)
from interfacea.io.atomrecord import AtomRecord
from interfacea.chemistry.elements import mapping as elem_map

logging.getLogger(__name__).addHandler(logging.NullHandler())


class PDBReader(object):
    """Yet another PDB Parser.

    This variation of a PDB format parser does not read lines other than
    coordinate information (ATOM, HETATM, MODEL, ENDMODEL). Like pretty
    much 99.9% of other PDB parsers out there.

    Args:
        filepath (Path): path to the PDB file to be read.
        permissive (bool, optional): ignores mal-formatted records in the PDB
            file if set to True. Defaults to False.
    Raises:
        FileNotFoundError: if the input file cannot be read.
        PDBReaderError: if a record line cannot be parsed correctly.
        PDBReaderWarning: if PDBReaderError is raised on an ATOM/HETATM line
            but permissive it set to True.
    """

    def __init__(self, filepath, permissive=False):

        self._line_parsers = {
            'ATOM': self._parse_atom,
            'HETATM': self._parse_atom,
            'MODEL': self._parse_model,
            'ENDMDL': self._parse_endmdl,
        }

        self._data = []

        self.path = filepath
        self.permissive = permissive
        self.read_file()
        self.check_models()

    @property
    def data(self):
        return self._data

    def read_file(self):
        """Parses the contents of a PDB-formatted file."""

        self._inmodel = False  # flag
        self._serial = 0  # we do not read serials from the file.
        self._model_no = 0  # we do not read model idxs from the file
        with self.path.open('r') as pdbfile:
            for lineno, line in enumerate(pdbfile, start=1):
                self.lineno, self.line = lineno, line
                self._parse_line()

        logging.debug(f"Read {self._serial + 1} atoms from file: {self.path}")

    def _parse_line(self):
        """Delegates parsing to appropriate methods"""

        self.rectype = self.line[:6].strip()
        try:
            p = self._line_parsers[self.rectype]
        except Exception:
            pass  # unsupported record
        else:
            p()  # parse line

    def _parse_model(self):
        """Parses MODEL records"""

        if self._inmodel:
            emsg = f"Missing ENDMDL record before line {self.lineno}"
            raise PDBReaderError(emsg) from None

        self._inmodel = True

    def _parse_endmdl(self):
        """Parses ENDMDL records"""
        if not self._inmodel:
            emsg = f"ENDMDL record outside of MODEL on line {self.lineno}"
            raise PDBReaderError(emsg) from None

        self._inmodel = False
        self._model_no += 1  # increment
        self._serial = 0  # reset

    def _parse_atom(self):
        """Parses ATOM/HETATM records"""

        try:
            line = self.line
            atom = AtomRecord(
                self._serial,
                self._model_no,
                self.rectype,
                line[12:16].strip(),
                line[16].strip(),
                line[17:20],
                line[21],
                int(line[22:26]),
                line[26].strip(),
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
                float(line[54:60]),
                float(line[60:66]),
                line[72:76].strip(),
                # line[76:78].strip(),  # auto-guess element
            )
        except Exception as err:
            emsg = f"Could not parse atom on line {self.lineno}"
            if self.permissive:
                emsg += " Ignoring."
                warnings.warn(emsg, PDBReaderWarning)
            else:
                raise PDBReaderError(emsg) from err
        else:
            atom.element = self._guess_element(line[12:16])
            self._data.append(atom)
            self._serial += 1

    def _guess_element(self, name_string):
        """Tries to infer the atomic element from the atom name.

        Adapted from Biopython's Atom.py code.

        Args:
            name_string (str): the full 4-char string defined in ATOM records.

        Returns:
            Element object.
        """

        stripped_name = name_string.strip()

        if name_string[0].isalpha() and not name_string[2:].isdigit():
            putative_elem = stripped_name
        elif stripped_name[0].isdigit():  # e.g. 1HE2
            putative_elem = stripped_name[1]
        else:
            putative_elem = stripped_name[0]

        element = elem_map.get(putative_elem.capitalize())
        if element.symbol is None:
            warnings.warn(
                f"Could not assign element from atom name: '{name_string}'",
                PDBReaderWarning
            )

        return element

    def check_models(self):
        """Asserts number of atoms is the same in all models.

        Raises:
            PDBReaderError: if the check fails.
        """

        a_p_m = collections.defaultdict(lambda: 0)  # atoms per model
        for r in self.data:
            a_p_m[r.model] += 1

        if len(set(a_p_m.values())) != 1:
            emsg = 'Models have different number of atoms'
            raise PDBReaderError(emsg)
