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
Unit tests for the utils module functions.
"""

from io import TextIOWrapper
import pathlib
import urllib.request as request

import pytest

from interfacea import utils

datadir = pathlib.Path('tests') / 'data'


def test_validate_existing_path():
    """Returns an existing path successfully"""

    p = datadir / '1brs.pdb'
    p_as_str = str(p)
    assert utils.validate_path(p_as_str).samefile(p)


def test_validate_nonexisting_path():
    """Raises error on non-existing path"""

    p = datadir / 'banana.pdb'
    p_as_str = str(p)
    with pytest.raises(FileNotFoundError, match='File not found:'):
        _ = utils.validate_path(p_as_str)


@pytest.mark.parametrize(
    ('code', 'fmt'),
    (
        ['1brs', 'pdb'],
    ),
)
@pytest.mark.usefixtures('has_internet_connectivity')
def test_fetch_rcsb_pdb(code, fmt):
    """Successfully loads a compressed PDB file from RCSB"""

    s = utils.fetch_rcsb_pdb(code, fmt=fmt)
    assert isinstance(s, TextIOWrapper)

    data = f'https://files.rcsb.org/download/{code}.pdb'  # uncompressed
    pdbtxt = request.urlopen(data)

    assert s.read() == pdbtxt.read().decode('utf-8')


@pytest.mark.usefixtures('has_internet_connectivity')
def test_fetch_rcsb_pdb_missing():
    """Raises error on RCSB 404 response"""

    code = 'arandomstringthatisnotapdbfilecodeforsure'

    with pytest.raises(ValueError, match='RCSB PDB entry'):
        _ = utils.fetch_rcsb_pdb(code)


@pytest.mark.usefixtures('has_internet_connectivity')
def test_fetch_rcsb_pdb_badfmt():
    """Raises error on unsupported fmt keyword"""

    code = '1brs'
    fmt = 'xyz'

    with pytest.raises(ValueError, match='File format not supported'):
        _ = utils.fetch_rcsb_pdb(code, fmt=fmt)
