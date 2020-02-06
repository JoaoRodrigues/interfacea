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
Unit tests for io.atomrecord sub module.
"""

import pytest

import interfacea.io.atomrecord as _ar


# atomrecord.AtomRecord
@pytest.mark.parametrize(
    "attrs",
    (
        # Passes
        {
            'serial': 0,
            'model': 1,
            'rectype': 'ATOM',
            'name': 'CA',
            'altloc': '',
            'resname': 'ALA',
            'chain': 'A',
            'resid': 15,
            'icode': '',
            'x': 1.0,
            'y': 1.0,
            'z': 1.0,
            'occ': 1.0,
            'b': 0.0,
        },
        {
            'serial': 0,
            'model': 1,
            'rectype': 'ATOM',
            'name': 'CA',
            'altloc': '',
            'resname': 'ALA',
            'chain': 'A',
            'resid': 15,
            'icode': '',
            'x': 1.0,
            'y': 1.0,
            'z': 1.0,
            'occ': 1.0,
            'b': 0.0,
            'segid': 'AAAA'
        },
        # Fails
        {
            'model': 1,
            'rectype': 'ATOM',
            'name': 'CA',
            'altloc': '',
            'resname': 'ALA',
            'chain': 'A',
            'resid': 15,
            'icode': '',
            'x': 1.0,
            'y': 1.0,
            'z': 1.0,
            'occ': 1.0,
            'b': 0.0,
        },  # missing serial
        {
            'serial': 'abc',
            'model': 1,
            'rectype': 'ATOM',
            'name': 'CA',
            'altloc': '',
            'resname': 'ALA',
            'chain': 'A',
            'resid': 15,
            'icode': '',
            'x': 1.0,
            'y': 1.0,
            'z': 1.0,
            'occ': 1.0,
            'b': 0.0,
        },  # serial is wrong type
    )
)
def test_make_AtomRecord(attrs):
    """Successfully creates AtomRecord class"""
    for trial in attrs:
        if 'serial' in attrs and isinstance(attrs['serial'], int):
            ar = _ar.AtomRecord(*attrs.values())
            for key, value in attrs.items():
                assert getattr(ar, key) == value

            if 'segid' not in attrs:  # default?
                assert ar.segid == ''

        else:
            with pytest.raises(TypeError):
                ar = _ar.AtomRecord(*attrs.values())
