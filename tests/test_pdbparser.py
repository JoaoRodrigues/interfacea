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
Unit tests for PDBParser.
"""

import pathlib
from urllib.request import urlopen

import pytest

from interfacea.io.pdb import (
    PDBParser,
    PDBParserError,
    # PDBParserWarning
)


datadir = pathlib.Path('tests') / 'data'


def test_no_arguments():
    """Throws error when creating PDBParser without arguments"""

    with pytest.raises(PDBParserError):
        _ = PDBParser()


def test_wrong_arguments():
    """Throws error when creating PDBParser with wrong arguments"""

    with pytest.raises(PDBParserError):
        _ = PDBParser(file='something', pdbid='else')


@pytest.mark.parametrize(
    'datafile',
    (
        datadir / 'default.pdb',
        datadir / 'multimodel.pdb'
    )
)
def test_load_local_files(datafile):
    """Create a new PDBParser from a file on disk"""

    p = PDBParser(datafile)
    assert p.source.read() == datafile.open('r').read()


@pytest.mark.usefixtures('has_internet_connectivity')
def test_load_remote_file():
    """Create a new PDBParser from a file on disk"""

    url = "https://files.rcsb.org/download/1ctf.pdb"
    p = PDBParser(pdbid='1ctf')
    assert p.source == urlopen(url).read().splitlines()


@pytest.mark.parametrize(
    'datafile',
    (
        datadir / 'default.pdb',
        datadir / 'multimodel.pdb'
    )
)
def test_parse(datafile):
    """Successfully parses test PDB files."""

    p = PDBParser(datafile)
    p.parse()

    # Assert some basic things manually.
    with datafile.open('r') as handle:
        natoms, nmodels = 0, 0
        for line in handle:
            if line.startswith('ATOM'):
                natoms += 1
            elif line.startswith('MODEL'):
                nmodels += 1

        if not nmodels and natoms:
            nmodels += 1

    assert len(p.atoms) == natoms // nmodels
    assert len(p.coords) > 0
    assert len(p.coords[0]) == nmodels


@pytest.mark.parametrize(
    ('datafile', 'emsg'),
    (
        (datadir / 'badatom.pdb', 'line 3'),
        (datadir / 'badmodel.pdb', 'Unexpected MODEL record'),
        (datadir / 'badendmdl.pdb', 'Unexpected ENDMDL record'),
        (datadir / 'duplicatedatom.pdb', 'is duplicated on line'),

    )
)
def test_parse_errors(datafile, emsg):
    """Successfully parses test PDB files."""

    p = PDBParser(datafile)
    with pytest.raises(PDBParserError, match=emsg):
        p.parse()
