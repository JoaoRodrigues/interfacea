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

import numpy as np

import interfacea.chemistry.elements as elements
from interfacea.core.atom import Atom
from interfacea.core.topology import Topology
from interfacea.core.structure import Structure
from interfacea.exceptions import InterfaceaError

from .helpers import guess_atom_element

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
        filepath = pathlib.Path(handle.name)

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
        topology : Topology
            topology containing atoms parsed from the PDB file.
        coordinates : list of lists
            list of x,y,z coordinates for each model in the parsed file.

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

        # Are we debug logging?
        self._print_debug = logger.isEnabledFor(logging.DEBUG)

        self._parse(handle)

    def _parse(self, handle):
        """Parse the file, line by line."""

        _recparsers = {
            "ATOM  ": self._parse_ATOM,
            "HETATM": self._parse_HETATM,
            "MODEL ": self._parse_MODEL,
            "ENDMDL": self._parse_ENDMDL,
            "END   ": self._terminate,
        }

        for ln, line in enumerate(handle, start=1):
            self.line = line
            self.ln = ln

            rectype = line[:6]

            try:
                _recparsers[rectype]()
            except KeyError:
                logger.debug(f"Ignoring record '{rectype}' on line {ln}")
                continue
            except StopIteration:
                break  # stop parsing file on demand
            except Exception as err:  # catch-all
                raise PDBParserError(f"Could not parse line {ln}") from err

        # Set coordinate array properly
        # in case we never found a MODEL statement
        if self.xyz and len(self.xyz[0]) == 3:  # only one model
            coordset = self.xyz[:]
            self.xyz = [coordset]

        self.xyz = np.array(self.xyz)

    # record parsers
    def _terminate(self):
        """Send stop signal to main parsing loop."""
        raise StopIteration

    def _parse_MODEL(self):
        """Reset the atomset, set the model flag, and add coordset."""

        if self._print_debug:
            logger.debug(f"Parsing MODEL statement at line {self.ln}")

        if self._in_model:
            raise ValueError(f"Unexpected MODEL record at line {self.ln}")

        self._atomset.clear()
        self._in_model = True

        self.xyz.append([])
        self._curxyz = self.xyz[-1]

    def _parse_ENDMDL(self):
        """Reset the model flag."""

        if self._print_debug:
            logger.debug(f"Parsing ENDMDL statement at line {self.ln}")

        if not self._in_model:
            raise PDBParserError(f"Unexpected ENDMDL record at line {self.ln}")

        self._in_model = False
        self._nrmodels += 1
        if self.nmodels is not None and self._nrmodels >= self.nmodels:
            self._terminate()

        # check all models have the same number of atoms
        if len({len(xyz) for xyz in self.xyz}) != 1:
            raise PDBParserError(
                f"Model '{self._nrmodels}' has a different number of atoms."
            )

    def _parse_coordinates(self, hetatm=False):
        """Parse an ATOM/HETATM record."""

        _guess_element = guess_atom_element

        if self._print_debug:
            logger.debug(f"Parsing coordinate data at line {self.ln}")

        line = self.line

        if not self._nrmodels:  # only parse topology info from first model
            fullname = line[12:16]
            try:
                elem = elements.ElementMapper[line[76:78].strip()]
            except KeyError:
                elem = _guess_element(fullname)
                if self._print_debug:
                    logger.debug(
                        f"Auto-assigned element {elem.symbol} for atom {fullname}"
                    )

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
