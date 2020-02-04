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

import pathlib

from interfacea.core import structure
from interfacea.io.base import Reader

__all__ = ['read']


# High-level IO
def read(filepath, **kwargs):
    """High-level method to create Structures from data files.

    Reader classes are picked based on the file extension, e.g. a .pdb file
    will be parsed by the io.PDBReader class. The mapping between extensions
    and readers is defined in io. Refer to that file and to each of the reader
    classes for more information on their arguments and options.

    Args:
        filepath (str): path to the file to be read.

    Returns:
        A Structure object containing atom coordinates and metadata.
    """

    # Validate Path
    try:
        path = pathlib.Path(filepath).resolve(strict=True)
    except FileNotFoundError as err:
        emsg = f"File not found or not readable: {filepath}"
        raise FileNotFoundError(emsg) from err
    except Exception as err:
        emsg = f"Unexpected error when parsing file path: {filepath}"
        raise IOError(emsg) from err

    # Parse Data
    r = Reader(path, **kwargs)
    name = kwargs.get('name', path.name)
    return structure.Structure.build(name, r.data, **kwargs)
