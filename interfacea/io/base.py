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

from interfacea.exceptions import InterfaceaError


class ParserError(InterfaceaError):
    pass


class StructureBuilderError(InterfaceaError):
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


class BaseStructureBuilder(metaclass=abc.ABCMeta):
    """Specifies an abstract interface for building Structure objects.

    Structure builders are essentially validators of the parsed data. They must
    implement a _build method that returns a Structure object.

    Methods
    -------
        build()
            Builds and returns a complete Structure.
    """

    @abc.abstractmethod
    def _build(self):
        """Builds an instance of Structure."""
        pass

    def build(self):
        """Returns a complete instance of a Structure"""
        raise NotImplementedError('Do this yourself!')
