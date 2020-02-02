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
Unit tests for io module classes.
"""

import pathlib

import pytest

import interfacea.io as _io
import interfacea.exceptions as e

@pytest.mark.parametrize(
    "ifile, altloc_flag",
    [
        (str(pathlib.Path("tests/data/pdb/default.pdb")), True),
        # (str(pathlib.Path("tests/data/pdb/default.pdb")), False),
    ]
)
def test_io_read(ifile, altloc_flag):
    """Successfully runs top-level read function"""

    _io.read(ifile, discard_altloc=altloc_flag)


@pytest.mark.parametrize(
    "ifile, error",
    [
        (str(pathlib.Path("tests/data/pdb/foo.bar")), FileNotFoundError),
        (1, Exception)
    ]
)
def test_io_read_fail(ifile, error):
    """Top-level read function fails on missing file"""

    with pytest.raises(error):
        _io.read(ifile)


class TestReader:
    """Tests for Reader class"""

    def setup(self):
        """Pre-test setup."""

        rootdir = pathlib.Path(".")  # rootdir
        self.datadir = rootdir / "tests" / "data"

    def test_read_known(self):
        """Successfully loads PDBReader"""

        _ = _io.base.Reader(self.datadir / "pdb" / "default.pdb")

    def test_unknown_extension(self):
        """Fails when extension is unknown/unsupported"""
        with pytest.raises(IOError):
            _ = _io.base.Reader(self.datadir / "pdb" / "foo.bar")


class TestPDBReader:
    """Tests for PDBReader class."""

    def setup(self):
        """Pre-test setup."""

        rootdir = pathlib.Path(".")  # rootdir
        self.datadir = rootdir / "tests" / "data" / "pdb"

    @pytest.mark.parametrize(
        "ifile, natoms, nmodels",
        [("default.pdb", 106, 1), ("multimodel.pdb", 212, 2)]
    )
    def test_read_pdb(self, ifile, natoms, nmodels):
        """Successfully load and parse a PDB file."""

        r = _io.pdb.PDBReader(self.datadir / ifile)
        atom_records = r.data
        assert len(atom_records) == natoms
        assert len({r.model for r in atom_records}) == nmodels

    @pytest.mark.parametrize(
        "ifile, errmessage",
        [
            ("nomodel.pdb", "ENDMDL record outside of MODEL on line 4"),
            ("noendmdl.pdb", "Missing ENDMDL record before line 4"),
            ("badatom.pdb", "Could not parse atom on line 3")
        ]
    )
    def test_read_pdb_fail(self, ifile, errmessage):
        """Fail when parsing a bad PDB (PDBFormatError)."""

        with pytest.raises(e.PDBFormatError) as excinfo:
            _io.pdb.PDBReader(self.datadir / ifile)

        assert excinfo.value.message == errmessage

    def test_permissive(self):
        """Ignore badly formatted lines when permissive=True"""

        with pytest.warns(e.PDBFormatWarning):
            _io.pdb.PDBReader(
                self.datadir / "badatom.pdb",
                permissive=True
            )
