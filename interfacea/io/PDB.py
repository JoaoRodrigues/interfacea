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

import logging
import warnings

from . import AtomRecord
from ..exceptions import (
    PDBFormatError,
    PDBFormatWarning
)

logging.getLogger(__name__).addHandler(logging.NullHandler())


class PDBReader(object):
    """Yet another PDB Parser.

    This variation of a PDB format parser does not read lines other than
    coordinate information (ATOM, HETATM, MODEL, ENDMODEL) and stops reading
    when an 'END' record is found. Like pretty much 99.9% of other PDB parsers
    out there.

    Args:
        filepath (str): path to the PDB file to be read.
        permissive (bool, optional): ignores mal-formatted records in the PDB
            file if set to True. Defaults to False.
    Raises:
        FileNotFoundError: if the input file cannot be read.
        PDBFormatError: if a record line cannot be parsed correctly.
        PDBFormatWarning: if PDBFormatError is raised on an ATOM/HETATM line
            but permissive it set to True.
    """

    def __init__(self, filepath, permissive=False):

        self.path = filepath
        self.permissive = permissive
        self._read_file()

    def read_file(self):
        """Parses the contents of a PDB-formatted file."""

        atom_list = []

        model = None
        serial = 0  # we do not read serials from the file.
        with self.path.open('rU') as pdbfile:
            for lineno, line in pdbfile:
                if line.startswith('MODEL'):  # we got ourselves an ensemble!
                    if model is not None:
                        emsg = f"PDB file is missing ENDMDL records before line {lineno}"
                        raise PDBFormatError(emsg) from None
                    try:
                        model = int(line[10:14])
                    except Exception:
                        emsg = f"Could not parse MODEL record on line {lineno}"
                        raise PDBFormatError(emsg) from None
                elif line.startswith(('ATOM', 'HETATM')):
                    if model is None:
                        model = 0
                    try:
                        atom = AtomRecord(
                            serial,
                            model,
                            line[:6].strip(),
                            line[12:16].strip(),
                            line[16].strip(),
                            line[17:20],
                            line[21],
                            int(line[22:26]),
                            line[26],
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                            float(line[54:60]),
                            float(line[60:66]),
                            line[72:76],
                            line[76:78],
                        )
                    except Exception as err:
                        emsg = f"Could not parse line no. {lineno}."
                        if self.permissive:
                            emsg += " Ignoring."
                            warnings.warn(emsg, PDBFormatWarning)
                        else:
                            raise PDBFormatError(emsg) from err
                    else:
                        atom_list.append(atom)
                        serial += 1

        logging.debug(f"Parsed {serial + 1} records from file: {self.filepath}")
        return atom_list
