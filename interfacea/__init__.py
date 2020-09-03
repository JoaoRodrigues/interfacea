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

from .config import set_log_level  # noqa: F401
from .constants import RANDOM_SEED  # noqa: F401
from .version import __version__  # noqa: F401

# Setup logger
# This is the parent logger since the library is supposed to be loaded from here.
# Hence, configs set here apply to all module-level loggers
logging.getLogger(__name__).addHandler(logging.NullHandler())
logging.getLogger(__name__).setLevel(logging.WARNING)
