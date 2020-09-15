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
Module containing classes to represent 3D molecular objects.

Atomic structures in interfacea are represented internally by the Structure
class, which stores both coordinate data and metadata. The metadata, e.g. atom
and residue names, is stored in individual Atom objects.
"""

import copy
import logging

import numpy as np

from interfacea.exceptions import InterfaceaError

# Setup module logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class StructureError(InterfaceaError):
    """Subclass for exceptions in this module."""

    pass


class Structure(object):
    """Represent 3D molecules as collections of atoms.

    Parameters
    ----------
        topology : Topology
            topology object describing the molecular structure.
        xyzdata : np.ndarray
            array with numerical coordinates of shape MxNx3, where M is the
            number of models (or frames) and N is the total number of atoms
            including altlocs.
        name : str, optional
            string that identifies the structure. Default is "Unnamed".

    Attributes
    ----------
        atoms : list
            ordered list of all atoms in the corresponding Topology.
        coords : np.ndarray
            coordinate array of the selected model.
        raw_coords : np.ndarray
            raw coordinate array including data for all models.
        nmodels : int
            number of models or frames in the structure.
    """

    def __init__(self, topology, xyzdata, name="Unnamed"):
        """Create a new instance of the class."""

        self.name = name

        self.topology = self._check_topology(topology)
        self._coords = self._check_coordinates(xyzdata)

        if self._coords.shape[1] != self.topology.nlocs:
            raise StructureError(
                "Topology size does not match coordinates: "
                f"{self._coords.shape[1]} != {self.topology.nlocs}"
            )

        self._model = 0
        self._kdtree = None

    # Internal dunder methods
    def __str__(self):
        """Return a string representation of the Structure."""
        return f"<Structure name='{self.name}' natoms={self.topology.natoms}>"

    def __iter__(self):
        """Return a generator to iterate over the children atoms."""
        yield from self.topology.atoms

    def __len__(self):
        """Return number of atoms."""
        return self.topology.natoms

    def __eq__(self, other):
        """Compare two structures by their topology and coordinates."""
        return self.topology == other.topology and \
            np.allclose(self._coords, other._coords)

    def __copy__(self):
        """Make shallow copy of structure."""
        return Structure(self.topology, self._coords)

    def __deepcopy__(self, memo):
        """Make deep copy of structure."""
        return Structure(copy.deepcopy(self.topology), copy.deepcopy(self._coords))

    # Private Methods
    def _check_topology(self, topology):
        """Check if the topology is valid."""

        # Check all atoms have unique indexes in the topology.
        idxset = {a.index for a in topology.unpack_atoms()}
        if len(idxset) != topology.nlocs:
            raise StructureError("Atoms in Topology must have unique indexes.")

        return topology

    def _check_coordinates(self, xyzdata):
        """Check if the coordinate array is valid.

        Arguments
        ---------
            xyzdata : list of lists or np.ndarray of floats
                array of shape MxNx3, where M is the number of models in the
                Structure and N is the number of atoms. The number of atoms
                _must_ match those in the corresponding Topology.
        """

        if not isinstance(xyzdata, np.ndarray):
            try:
                xyzdata = np.array(xyzdata, dtype=np.float)
            except Exception as err:
                raise StructureError("Cannot coerce coordinates to np.array") from err

        if xyzdata.ndim != 3:
            raise StructureError(
                f"Coordinate array has wrong dimensions: {xyzdata.ndim} != 3"
            )

        if xyzdata.shape[2] != 3:
            raise StructureError(
                f"Wrong coordinate array shape {xyzdata.shape} != (x, x, 3)"
            )

        return xyzdata

    # Public Methods
    @property
    def model(self):
        """Return the index of the active model."""
        return self._model

    @model.setter
    def model(self, index):
        """Define which model to pick data from in multi-model structures.

        Arguments
        ---------
            index : int
                (0-based) index of a model to select.
        """

        if not 0 <= index < self.nmodels:
            emsg = f"Model index must be between 0 and {self.nmodels - 1}"
            raise StructureError(emsg)

        self._model = index
        self._kdtree = None  # delete cached kdtree
        logging.debug(f"Set active model to: {index}")

    @property
    def nmodels(self):
        """Return the number of models in the structure."""
        return self._coords.shape[0]

    @property
    def coords(self):
        """Return a view of the coordinate array for the selected model."""
        return self._coords[self._model, :]

    @property
    def raw_coords(self):
        """Return the entire coordinate array (all models)."""
        return self._coords
