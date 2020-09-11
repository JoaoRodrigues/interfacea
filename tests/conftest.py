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

"""General Settings for Unit Tests."""

import pathlib
import socket

import pytest


# Connectivity Fixture
@pytest.fixture(scope="session")
def has_internet_connectivity():
    """Return True if google can be reached on port 80."""
    try:
        socket.create_connection(("www.google.com", 80))
        return
    except OSError:
        pass

    pytest.skip("no internet connection")


# Define data folder relative to this file's location
@pytest.fixture(scope="session")
def datadir():
    """Define data folder."""
    this = pathlib.Path(__file__).resolve().parent
    return this / "data"
