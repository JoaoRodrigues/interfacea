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
Module containing classes to parse structural data files.
"""

from interfacea.exceptions import InterfaceaError
import interfacea.utils as ia_utils

# Bring parsers to this namespace for convenience
from .pdb import PDBParser

# Bring builders to this namespace for convenience
from .aabuilder import AllAtomStructureBuilder


_parsers = {
    'pdb': PDBParser,
    'ent': PDBParser,
}


# Module-level functions
def read(filepath, name=None, parser=None, builder=None):
    """High-level method to create Structure objects from files.

    The parser and builder arguments define how the file is parsed and how the
    Structure is built. Refer to the appropriate modules for more details. By
    default, the parser is picked based on the input file extension and the
    AllAtomStructureBuilder is used to build the structure.

    Arguments
    ---------
        filepath : str
            path to the file to be read.

        name : str, optional
            identifier for the Structure. If None, the input file name will be
            used as a name.

        parser : BaseParser, optional
            format-specific parser to use when reading the input file.

        builder : BaseStructureBuilder, optional
            factory builder class to assemble the structure.

    Returns
    -------
        a Structure object containing atom coordinates and metadata.
    """

    fpath = ia_utils.validate_path(filepath)
    if parser is None:
        try:
            fext = fpath.suffix.lstrip('.')
            parser = _parsers[fext](fpath.open('rt'))
        except KeyError:
            emsg = (
                f'Could not auto-assign a parser to file "{filepath}"'
                f'\nSupported file types: {",".join(_parsers)}'
            )
            raise InterfaceaError(emsg)

    if builder is None:
        builder = AllAtomStructureBuilder

    if name is None:
        name = str(fpath.stem)

    b = builder(parser)
    structure = b.build()
    structure.name = name

    return structure


def fetch(pdbcode, data_fmt='pdb', name=None, parser=None, builder=None):
    """High-level method to create Structure objects from remote data.

    Arguments
    ---------
        pdbcode : str
            4-character RCSB entry code.

        data_fmt : str, optional
            Format in which to download the data from RCSB.

        name : str, optional
            identifier for the Structure. If None, the pdb code will be used as
            a name.

        parser : BaseParser, optional
            format-specific parser to use when reading the data.

        builder : BaseStructureBuilder, optional
            factory builder class to assemble the structure.

    Returns
    -------
        a Structure object containing atom coordinates and metadata.
    """

    gzfile = ia_utils.fetch_rcsb_pdb(pdbcode, data_fmt)
    if parser is None:
        parser = _parsers[data_fmt](gzfile)

    if builder is None:
        builder = AllAtomStructureBuilder

    if name is None:
        name = pdbcode

    b = builder(parser)
    structure = b.build()
    structure.name = name

    return structure
