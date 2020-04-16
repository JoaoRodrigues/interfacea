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
Module containing assorted utility functions shared across modules.
"""

import gzip
import io
import logging
import pathlib
from urllib import (
    error as url_error,
    request as url_request
)

logging.getLogger(__name__).addHandler(logging.NullHandler())


def validate_path(path):
    """Asserts that the path is valid and exists.

    Arguments
    ---------
        path : str
            path to a file on disk, as a string.

    Returns
    -------
        the full path as a pathlib.Path object.
    """

    p = pathlib.Path(path)
    try:
        return p.resolve(strict=True)
    except FileNotFoundError:
        emsg = f"File not found: {path}"
        raise FileNotFoundError(emsg) from None


def fetch_rcsb_pdb(pdbid, fmt='pdb'):
    """Fetches a structure from the RCSB PDB database.

    Always downloads the full assymetric unit.

    Arguments
    ---------
        pdbid : str
            4-character RCSB PDB code as a string.

        fmt : str, optional
            file type to download: 'pdb' (default) or 'cif'.

    Returns
    -------
        structure data as a file-like object
    """

    _supported_fmt = set(('pdb', 'cif'))
    if fmt not in _supported_fmt:
        raise ValueError(
            f'File format not supported: {fmt}'
            f'\nChoose from: {",".join(_supported_fmt)}'
        )

    rooturl = "https://files.rcsb.org/download/"
    fullurl = rooturl + pdbid + '.' + fmt + '.gz'  # get compressed by default

    logging.debug(f'Reading remote file from: {fullurl}')

    try:

        with url_request.urlopen(fullurl) as response:
            gz_data = io.BytesIO(response.read())
            gzfile = gzip.open(gz_data, mode='rt', encoding='utf-8')
            return gzfile

    except url_error.URLError as err:
        if err.code == 404:
            raise ValueError(f'RCSB PDB entry {pdbid} does not exist.')
        raise  # re-raise
