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


class StructureBuilderError(BaseInterfaceaException):
    pass


class BaseParser(metaclass=abc.ABCMeta):
    """Specifies an abstract interface for parsing structure files.

    Attributes
    ----------
        atoms : list of Atom
            sequence of Atom objects representing the information contained in
            the file.

        coords : list of list of float
            nested structure of shape MxNx3 representing all cartesian
            coordinates read from the file, where M is the number of models and
            N is the number of atoms.
    """

    @abc.abstractmethod
    def parse(self):
        """Parses the content of the file and returns the atoms/coords"""
        pass

    @property
    def atoms(self):
        return self._atoms

    @property
    def coords(self):
        return self._coords


# class BaseStructureBuilder(metaclass=abc.ABCMeta):
#     """Specifies an abstract interface for building Structure objects.
#     """

#     def __init__(self, atoms, coords):
#         pass

#     @abc.abstractmethod
#     def build(self):
#         """Returns an instance of Structure."""
#         return self.structure


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
