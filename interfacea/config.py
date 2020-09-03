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

"""Module containing assorted utility functions shared across modules."""

import logging


# Functions
def set_log_level(level="verbose"):
    """Define the verbosity of the logging messages.

    Parameters
    ----------
        level : str
            Verbosity of logging. Choose from minimal, verbose, or debug. The
            default setting is not to setup any logging. You have to call this
            function explicitly to see messages other than warnings or errors.

    Raises
    ------
        ValueError
            The logging level you provided is not supported.
    """

    # Logging levels
    _level_dict = {"minimal": 30, "verbose": 20, "debug": 10}

    log_level = _level_dict.get(level)
    if log_level is None:
        raise ValueError(f"Unsupported log level: {level}")

    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        fmt="[%(asctime)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.handlers = []  # clear handler list

    root_logger.addHandler(handler)
    root_logger.setLevel(log_level)

    logging.critical(f"Logging enabled (level={level})")  # always show
