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
Exceptions raised by modules and classes across interfacea.
"""

# Warnings


class PDBReaderWarning(Warning):
    pass


# Exceptions


class InterfaceBaseException(Exception):
    """Base class for exceptions in this module"""

    def __init__(self, message):
        self.message = message


class DisorderedAtomError(InterfaceBaseException):
    pass


class PDBReaderError(InterfaceBaseException):
    pass


class FunctionalGroupError(InterfaceBaseException):
    pass


class StructureBuilderError(InterfaceBaseException):
    pass
