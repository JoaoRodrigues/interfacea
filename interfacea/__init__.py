#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018 João Pedro Rodrigues
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
The interfacea package contains tools to analyze features of protein interfaces.
"""

import logging
import random

from .io import read  # noqa: F401
from .version import __version__  # noqa: F401

# Setup logger
# This is the parent logger since the library is supposed
# to be loaded from here. Hence, configs set here apply to
# all module-level loggers
logging.getLogger(__name__).addHandler(logging.NullHandler())
logging.getLogger(__name__).setLevel(logging.WARNING)

# Global Constants
RANDOM_SEED = random.randint(0, 1000)  # user can set it manually later


# Global Methods
def set_log_level(level='verbose'):
    """Controls verbosity of the logs.

    Useful for interactive/debugging applications.

    Args:
        level (str): verbosity/type of logging. Choose from:
            - 'silent': disables logging.
            - 'minimal': only warnings and other critical messages
            - 'verbose': informative/descriptive messages (default).
            - 'debug': very verbose internal/diagnostic messages.

    Raises:
        ValueError: logging level is not supported.
    """

    # Logging levels
    _level_dict = {
        'minimal': logging.WARNING,
        'verbose': logging.INFO,
        'debug': logging.DEBUG
    }

    log_level = _level_dict.get(level)
    if log_level is None:
        emsg = f"Unsupported log level: {level}"
        raise ValueError(emsg)

    handler = logging.StreamHandler()
    formatter = logging.Formatter(fmt='[%(asctime)s] %(message)s',
                                  datefmt='%H:%M:%S')
    handler.setFormatter(formatter)

    # Override root logger - assume we only call this function interactively.
    root_logger = logging.getLogger()
    root_logger.handlers = []  # clear handler list
    root_logger.addHandler(handler)

    root_logger.setLevel(log_level)
    logging.critical(f"Logging enabled (level={level})")  # always show


# Randomness
def set_random_seed(seed=917):
    """Defines the numerical random seed for reproducibility.

    Setting this seed does not ensure *complete reproducibility*. Some methods
    are not deterministic across different hardware configurations even with
    the same random seed.

    Args:
        seed (int): positive integer to use as a random seed.

    Raises:
        TypeError: seed number is not a valid positive number.
    """

    if not (isinstance(seed, int) and seed > 0):
        raise TypeError(f"Seed must be a positive integer: {seed}")

    global RANDOM_SEED
    RANDOM_SEED = seed
