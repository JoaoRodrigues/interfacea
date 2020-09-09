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

"""Module containing classes to build molecular structures from PDB files."""

import gzip
import logging
from io import TextIOBase
import pathlib
import warnings

import interfacea.chemistry.elements as elements
from interfacea.core.atom import Atom
from interfacea.core.topology import Topology
from interfacea.core.structure import Structure
from interfacea.exceptions import InterfaceaError

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class PDBParserError(InterfaceaError):
    """Exception specific to this submodule."""

    pass


class PDBParserWarning(UserWarning):
    """Warning class specific to this submodule."""

    pass


class PDBParser:
    """Yet another PDB Parser.

    This very original variation of a PDB format parser reads only coordinate
    records: ATOM, HETATM, MODEL, ENDMDL. Like 99.9% of all other PDB parsers
    out there.

    The parser is strict (quits on badly-formatted lines) and does very little
    validation of what it does parse successfully. It is the user"s job to make
    sure that the input data actually makes sense: ie. garbage in, garbage out.

    Arguments
    ---------
        handle: TextIOBase
            Open file or file-like object containing the PDB file data in text
            format.

    Attributes
    ----------
        structure : Structure
            Structure object representing the parsed PDB file.

    Raises
    ------
        PDBParserError
    """

    def __init__(self, handle):
        """Create a new instance of the class."""
        self._top = Topology()
        self._xyz = []
        self.structure = None

        self._atomset = set()  # to avoid duplicates in each model
        self._in_model = False  # flag set when reading multi-model files.
        self._nmodels = 0

        if not isinstance(handle, TextIOBase):
            raise PDBParserError(f"Wrong type for handle argument: {type(handle)}")

        self.parse(handle)

    @staticmethod
    def from_file(filepath):
        """Parse the atomic information in the PDB file into Atom objects.

        Arguments
        ---------
            filepath : pathlib.Path or str
                location of the input PDB file.

        Returns
        -------
            Structure object.
        """

        if isinstance(filepath, str):
            filepath = pathlib.Path(filepath)
        elif isinstance(filepath, pathlib.Path):
            pass
        else:
            raise TypeError(f"Unsupported input type: {type(filepath)}")

        filepath = filepath.resolve(strict=True)

        if "gz" in filepath.suffix:
            opener = gzip.open(filepath, "rt")
        else:
            opener = filepath.open("r")

        with opener as handle:
            reader = PDBParser(handle)

        return Structure(name=filepath.name, topology=reader._top, xyzdata=reader._xyz)

    def parse(self, handle):
        """Parse the file, line by line."""

        for ln, line in enumerate(handle, start=1):
            self.line = line
            self.ln = ln

            rectype = line[:6].strip()

            try:
                p = getattr(self, f"_parse_{rectype}")
            except AttributeError:
                logging.debug(f"Ignoring record '{rectype}' on line {ln}")
            else:
                try:
                    p()
                except Exception as err:  # catch-all
                    raise PDBParserError(f"Could not parse line {ln}") from err

        # Set coordinate array properly
        # in case we never found a MODEL statement
        if self._xyz and len(self._xyz[0]) == 3:  # only one model
            coordset = self._xyz[:]
            self._xyz = [coordset]

    # record parsers
    def _parse_MODEL(self):
        """Reset the atomset, set the model flag, and add coordset."""

        logging.debug(f"Parsing MODEL statement at line {self.ln}")

        if self._in_model:
            raise ValueError(f"Unexpected MODEL record at line {self.ln}")

        self._atomset.clear()
        self._in_model = True

        if self._xyz and len(self._xyz[0]) == 3:  # first model
            coordset = self._xyz[:]
            self._xyz = [coordset]

        self._xyz.append([])

    def _parse_ENDMDL(self):
        """Reset the model flag."""

        logging.debug(f"Parsing ENDMDL statement at line {self.ln}")

        if not self._in_model:
            raise PDBParserError(f"Unexpected ENDMDL record at line {self.ln}")

        self._in_model = False
        self._nmodels += 1

        # check all models have the same number of atoms
        atoms_in_models = {len(xyz) for xyz in self._xyz}
        assert (
            len(atoms_in_models) == 1
        ), f"Model '{self._nmodels} has a different number of coordinates."

    def _parse_coordinates(self, hetatm=False):
        """Parse an ATOM/HETATM record."""

        logging.debug(f"Parsing coordinate data at line {self.ln}")

        line = self.line

        if not self._nmodels:  # only parse topology info from first model
            fullname = line[12:16]
            try:
                elem = elements.ElementMapper[line[76:78].strip()]
            except KeyError:
                elem = self._guess_element(fullname)

            metadata = dict(
                hetatm=hetatm,
                altloc=line[16].strip(),
                resname=line[17:20],
                chain=line[21],
                resid=int(line[22:26]),
                icode=line[26].strip(),
                occ=float(line[54:60]),
                b=float(line[60:66]),
                segid=line[72:76].strip(),
                element=elem,
            )

            # build Atom object and add to Topology
            atom = Atom(name=fullname.strip(), **metadata)
            self._top.add_atom(atom)

        # Parse coordinates
        coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        self._xyz[-1].append(coords)

    def _parse_ATOM(self):
        """Parse an ATOM record."""
        self._parse_coordinates()

    def _parse_HETATM(self):
        """Parse a HETATM record."""
        self._parse_coordinates(hetatm=True)

    # Auxiliary Methods
    def _guess_element(self, fullname):
        """Infer the atomic element from the atom name.

        Adapted from Biopython"s Atom.py code.

        Arguments
        ---------
            fullname : str
                the full 4-char string defined in the PDB file.

        Returns
        -------
            Element object or elements.Unknown if ambiguous/unsuccessful.

        """

        stripped_name = fullname.strip()

        if fullname[0].isalpha() and not fullname[2:].isdigit():
            putative_elem = stripped_name
        elif stripped_name[0].isdigit():  # e.g. 1HE2
            putative_elem = stripped_name[1]
        else:
            putative_elem = stripped_name[0]

        element = elements.ElementMapper.get(putative_elem.capitalize())
        if element is None:
            warnings.warn(
                f"Could not assign element for atom {repr(fullname)} at line {self.ln}",
                PDBParserWarning,
            )
            element = element.Unknown

        logging.debug(f"Auto-assigned element {element.symbol} for atom {fullname}")
        return element
