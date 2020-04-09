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

import gzip
import logging
import pathlib
from urllib.request import urlopen
import warnings

from interfacea.chemistry.elements import ElementMapper, Unknown
from interfacea.core.atom import Atom
from interfacea.io.base import BaseParser, ParserError

logging.getLogger(__name__).addHandler(logging.NullHandler())


class PDBParserError(ParserError):
    pass


class PDBParserWarning(UserWarning):
    pass


class PDBParser(BaseParser):
    """Yet another PDB Parser.

    This very original variation of a PDB format parser reads only coordinate
    records: ATOM, HETATM, MODEL, ENDMDL. Like 99.9% of all other PDB parsers
    out there.

    The parser is strict (quits on badly-formatted lines) but does very little
    validation of what it does parse successfully. It is the user's job to make
    sure that the input data actually makes sense: ie. GIGO.

    Arguments
    ---------
        file : str, optional
            path to a PDB file stored locally.

        pdbid : str, optional
            four-character RCSB PDB id for the structure, e.g. 1brs

    Attributes
    ----------
        atoms : list of Atom
            list of Atom objects representing unique atoms present in the
            file. For multi-model files, only the atoms of the first model are
            kept.

        coords : list of list of float
            cartesian coordinates of every atom across all models of the PDB
            file. For multi-model files, coords has shape NxMx3, where N is the
            number of atoms and M is the number of models. For single structure
            files, M=1.

        source : object
            file-like object from where the data was read.

    Raises
    ------
        PDBParserError
    """

    def __init__(self, file=None, pdbid=None, stream=None):

        self.atoms = []
        self.coords = []

        self._atomdict = {}  # maps atom full id to serial to allocate coords
        self._atomset = set()  # to avoid duplicates in each model
        self._in_model = False  # flag set when reading multi-model files.

        if (file == pdbid is None) or (file is not None and pdbid is not None):
            raise PDBParserError('Provide either the file or pdbid argument')
        elif file is not None:
            self.source = self._validate_path(file)
        elif pdbid is not None:
            self.source = self._fetch_rcsb_pdb(pdbid)

    def _validate_path(self, path):
        """Asserts that the path is valid and exists.

        Arguments
        ---------
            path : str
                path to a file, as a string.
        """

        logging.debug(f'Reading file from disk: {path}')

        p = pathlib.Path(path)
        try:
            return p.resolve(strict=True).open('r')
        except FileNotFoundError:
            emsg = f"File not readable: {path}"
            raise ParserError(emsg)

    def _fetch_rcsb_pdb(self, pdbid):
        """Fetches a PDB file from the RCSB PDB database.

        Always downloads the full assymetric unit.
        """

        rooturl = "https://files.rcsb.org/download/"
        fullurl = rooturl + pdbid + '.pdb.gz'  # get compressed by default

        logging.debug(f'Reading remote file from: {fullurl}')

        response = urlopen(fullurl)
        return gzip.decompress(response.read()).splitlines()

    def parse(self):
        """Parses the atomic information in the PDB file into Atom objects."""

        # For each atom line, create an Atom object.
        # Put coordinates in _x, _y, _z attributes (lists)
        for n, line in enumerate(self.source, start=1):
            self._process_line(n, line)

        return (self.atoms, self.coords)

    def _process_line(self, n, line):
        """Parses an individual line of the PDB file.

        Arguments
        ---------
            n : int
                line number in the original file.
            line : str
                contents of the line, as a string.
        """

        self.n = n
        self.line = line

        rectype = line[:6].strip()
        try:
            p = getattr(self, f'_parse_{rectype}')
        except AttributeError:
            logging.debug(f'Ignoring record "{rectype}" on line {n}')
        else:
            try:
                p()
            except PDBParserError:
                raise  # re-throw intact
            except Exception as err:  # catch-all
                raise PDBParserError(
                    f'Could not parse line {self.n}'
                ) from err

    # record parsers
    def _parse_MODEL(self):
        """Resets the atomdict and sets the model flag"""

        logging.debug(f'Parsing MODEL statement at line {self.n}')

        if self._in_model:
            raise PDBParserError(f'Unexpected MODEL record at line {self.n}')

        self._atomset.clear()
        self._in_model = True

    def _parse_ENDMDL(self):
        """Resets the model flag"""

        logging.debug(f'Parsing ENDMDL statement at line {self.n}')

        if not self._in_model:
            raise PDBParserError(f'Unexpected ENDMDL record at line {self.n}')

        self._in_model = False

    def _parse_ATOM(self):
        """Parses an ATOM record."""

        logging.debug(f'Parsing ATOM statement at line {self.n}')

        line = self.line

        fullname = line[12:16]
        elem = line[76:78].strip()
        if not elem:
            elem = self._guess_element(fullname)

        metadata = dict(
            hetatm=False,
            altloc=line[16].strip(),
            resname=line[17:20],
            chain=line[21],
            resid=int(line[22:26]),
            icode=line[26].strip(),
            occ=float(line[54:60]),
            b=float(line[60:66]),
            segid=line[72:76].strip(),
            element=elem
        )

        # build Atom object
        atom = Atom(
            name=fullname.strip(),
            serial=len(self.atoms),
            **metadata
        )

        # Parse coordinates
        coords = [
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ]

        self._commit_atom(atom, coords)

    def _parse_HETATM(self):
        """Parses a HETATM record."""

        self._parse_ATOM()  # how lazy can we be?
        self.atoms[-1].hetatm = True

    # Auxiliary Methods
    def _commit_atom(self, atom, coords):
        """Validates the Atom and appends it to the parsed data if successful"""

        self._check_duplicate(atom)

        serial = self._atomdict.get(atom.full_id)
        if serial is None:
            logging.debug(f'Registering new atom in PDBParser: {atom.full_id}')
            serial = self._atomdict[atom.full_id] = atom.serial

            self.coords.append([])
            self.atoms.append(atom)

        self._atomset.add(atom.full_id)
        self.coords[serial].append(coords)

    def _check_duplicate(self, atom):
        """Ensures the atom has not been read yet."""

        if atom.full_id in self._atomset:
            raise PDBParserError(
                f'Atom {atom} is duplicated on line {self.n}'
            )

    def _guess_element(self, fullname):
        """Tries to infer the atomic element from the atom name.

        Adapted from Biopython's Atom.py code.

        Args:
            fullname (str): the full 4-char string defined in the PDB file.

        Returns:
            Element object.
        """

        stripped_name = fullname.strip()

        if fullname[0].isalpha() and not fullname[2:].isdigit():
            putative_elem = stripped_name
        elif stripped_name[0].isdigit():  # e.g. 1HE2
            putative_elem = stripped_name[1]
        else:
            putative_elem = stripped_name[0]

        element = ElementMapper[putative_elem.capitalize()]
        if element is Unknown:
            warnings.warn(
                f'Could not assign element for atom {repr(fullname)} '
                f'at line {self.n}',
                PDBParserWarning
            )

        logging.debug(
            f'Auto-assigned element {element.symbol} for atom {fullname}'
        )
        return element
