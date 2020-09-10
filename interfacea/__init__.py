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

"""Top-level module of the interfacea package to analyze protein interfaces."""

import logging

from .constants import RANDOM_SEED  # noqa
from .version import __version__  # noqa
from .io import read  # noqa

# Setup logger
# This is the parent logger since the library is supposed to be loaded from here.
# Hence, configs set here apply to all module-level loggers
logging.basicConfig(
    format="[%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.WARNING,
)


# Logging control helpers
def verbose():
    """Set loggers to verbose."""
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.critical("Logging set to verbose")


def default():
    """Set loggers to default."""
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    logging.critical("Logging set to default")


def debug():
    """Set loggers to debug."""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logging.critical("Logging set to debug")


def silent():
    """Set silent _most_ logging messages."""
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)
    logging.critical("Logging set to silent")
