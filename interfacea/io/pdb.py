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


def read_pdb(filepath, **kwargs):
    """Produce a Structure from a PDB file.

    Arguments
    ---------
        filepath : TextIOBase or str
            handle or path to the file to be parsed.

        nmodels : int, optional
            number of models to read from the file. Default, all.

    Returns
    -------
        a Structure object containing a topology and coordinates.
    """

    if isinstance(filepath, TextIOBase):
        handle = filepath
        filepath = pathlib.Path(handle)

    elif isinstance(filepath, str):
        filepath = pathlib.Path(filepath)

        _ = filepath.resolve(strict=True)  # to validate that path exists

        if filepath.suffix == ".gz":
            handle = gzip.open(filepath, "rt")
        else:
            handle = filepath.open("r")

    elif isinstance(filepath, pathlib.Path):
        _ = filepath.resolve(strict=True)  # to validate that path exists
        if filepath.suffix == ".gz":
            handle = gzip.open(filepath, "rt")
        else:
            handle = filepath.open("r")
    else:
        raise ValueError(f"Unknown argument type for _open_file: {type(filepath)}")

    with handle:
        reader = PDBParser(handle, **kwargs)

    return Structure(name=filepath.name, topology=reader.topology, xyzdata=reader.xyz)


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
        nmodels : int, optional
            number of models to read from the file. Default is None, meaning
            read all frames.

    Attributes
    ----------
        structure : Structure
            Structure object representing the parsed PDB file.

    Raises
    ------
        PDBParserError
    """

    def __init__(self, handle, nmodels=None):
        """Create a new instance of the class."""

        self.topology = Topology()
        self.xyz = []
        self._curxyz = self.xyz  # pointer

        self.structure = None
        self.nmodels = nmodels

        self._atomset = set()  # to avoid duplicates in each model
        self._in_model = False  # flag set when reading multi-model files.
        self._nrmodels = 0  # number of READ models.

        if not isinstance(handle, TextIOBase):
            raise PDBParserError(f"Wrong type for handle argument: {type(handle)}")

        self.parse(handle)

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
                except StopIteration:
                    break  # stop parsing file on demand
                except Exception as err:  # catch-all
                    raise PDBParserError(f"Could not parse line {ln}") from err

        # Set coordinate array properly
        # in case we never found a MODEL statement
        if self.xyz and len(self.xyz[0]) == 3:  # only one model
            coordset = self.xyz[:]
            self.xyz = [coordset]

    # record parsers
    def _parse_MODEL(self):
        """Reset the atomset, set the model flag, and add coordset."""

        logging.debug(f"Parsing MODEL statement at line {self.ln}")

        if self._in_model:
            raise ValueError(f"Unexpected MODEL record at line {self.ln}")

        self._atomset.clear()
        self._in_model = True

        if self.xyz and len(self.xyz[0]) == 3:  # first model
            coordset = self.xyz[:]
            self.xyz = [coordset]

        self.xyz.append([])
        self._curxyz = self.xyz[-1]

    def _parse_ENDMDL(self):
        """Reset the model flag."""

        logging.debug(f"Parsing ENDMDL statement at line {self.ln}")

        if not self._in_model:
            raise PDBParserError(f"Unexpected ENDMDL record at line {self.ln}")

        self._in_model = False
        self._nrmodels += 1
        if self.nmodels is not None and self._nrmodels >= self.nmodels:
            raise StopIteration  # stop parsing

        # check all models have the same number of atoms
        atoms_in_models = {len(xyz) for xyz in self.xyz}
        assert (
            len(atoms_in_models) == 1
        ), f"Model '{self._nrmodels} has a different number of coordinates."

    def _parse_coordinates(self, hetatm=False):
        """Parse an ATOM/HETATM record."""

        logging.debug(f"Parsing coordinate data at line {self.ln}")

        line = self.line

        if not self._nrmodels:  # only parse topology info from first model
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
                occupancy=float(line[54:60]),
                b=float(line[60:66]),
                segid=line[72:76].strip(),
                element=elem,
            )

            # build Atom object and add to Topology
            atom = Atom(name=fullname.strip(), **metadata)
            self.topology.add_atom(atom)

        # Parse coordinates
        coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        self._curxyz.append(coords)

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
