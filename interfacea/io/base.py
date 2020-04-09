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
Base class for file parsers.
"""

import abc

from interfacea.exceptions import BaseInterfaceaException


class ParserError(BaseInterfaceaException):
    pass


class BaseParser(metaclass=abc.ABCMeta):
    """Specifies an abstract interface for parsing structure files.
    """

    def __init__(self, filepath, **kwargs):

        self.fpath = self._validate_path(filepath)

    @abc.abstractmethod
    def parse(self):
        """Parses the content of the file and returns a list of objects"""
        pass


# class Reader(object):
#     """Delegates parsing of a file to an appropriate reader class."""

#     def __new__(cls, path, **kwargs):
#         try:
#             r = readers[path.suffix]
#         except KeyError:
#             emsg = f"File format not supported ({path.suffix})"
#             raise IOError(emsg) from None
#         else:
#             return r(path, kwargs)
