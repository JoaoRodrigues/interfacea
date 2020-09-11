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

"""Module containing classes to parse structural data files."""

import logging
import pathlib
import warnings

from interfacea.core.topology import Topology
from interfacea.exceptions import InterfaceaError, InterfaceaWarning
from .pdb import read_pdb

__all__ = ["read"]

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

_READERS = {
    ".pdb": read_pdb,
    ".pdb.gz": read_pdb,
    ".ent": read_pdb,
    ".ent.gz": read_pdb,
}

_INCLUDE_TOPOLOGY = {  # set of formats that include topology information
    ".pdb",
    ".pdb.gz",
    ".ent",
    ".ent.gz",
}


def _get_extension(filepath):
    """Return the extension from a string containing a file path."""

    if not isinstance(filepath, pathlib.Path):
        filepath = pathlib.Path(filepath)

    p = filepath.resolve(strict=True)

    fext = p.suffix
    if fext == ".gz":
        fext = "".join(p.suffixes[-2:])  # catch pdb.gz
    logger.info(f"Loading file '{str(p)}' with extension '{fext}'")
    return fext


def read(filepath, **kwargs):
    """High-level method to create Structure objects from files.

    Arguments
    ---------
        filepath : str
            path to the file to be parsed.

        topology : Topology or str, optional
            path to the file to be read as a topology or a
            pre-built Topology object. Required if the main input file does
            not define a topology (e.g. trajectories).

    Returns
    -------
        a Structure object containing a topology and coordinates.
    """

    fext = _get_extension(filepath)

    if kwargs.get("topology") is not None:
        if fext in _INCLUDE_TOPOLOGY:
            warnings.warn(
                "Ignoring provided topology as file includes topology information",
                InterfaceaWarning,
            )
            kwargs.pop("topology")
        else:
            if isinstance(kwargs["topology"], Topology):
                pass
            else:
                topext = _get_extension(kwargs["topology"])
                r = _READERS[topext](kwargs["topology"], nmodels=0)
                kwargs["topology"] = r.topology

    try:
        parser = _READERS[fext]
    except KeyError:
        emsg = (
            f"Could not auto-assign a parser to file '{filepath}' with extension '{fext}'"
            f"\nSupported file types: {','.join(_READERS)}"
        )
        raise InterfaceaError(emsg)
    else:
        return parser(filepath, **kwargs)
