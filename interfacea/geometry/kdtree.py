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

"""Wrapper code for KDTree C implementation."""

import numpy as np

from .kdtrees import KDTree as CKDTree


class KDTree:
    """Wrapper class around KDTree C code."""

    def __init__(self, structure, bucket_size=10):
        """Create a new KDTree instance from a Structure.

        Arguments
        ---------
            structure : Structure
                a molecular structure with topology and coordinates.
            bucket_size : int
                bucket size of the KDTree, or the maximum number of data points
                stored under in leaf node of the tree. Default is 10.
        """

        if bucket_size < 1:
            raise ValueError(f"Bucket size must be a positive integer: '{bucket_size}'")

        self.kdtree = CKDTree(structure.coords, bucket_size)
        self.structure = structure

    def find_neighbors(self, center, radius, sort=False, distances=False):
        """Return all Atom neighbors within radius of the center.

        Arguments
        ---------
            center: Atom or nd.array
                central point from which to find neighbors. Can be either an
                Atom object or a 3-item array. If the input is a DisorderedAtom,
                the coordinates of the selected child will be used as the center.
            radius : float
                distance threshold to consider an Atom neighbor of the center
                point. Same units as the structure coordinates.
            sort : bool
                sorts the neighbors in ascending order of distance to center,
                i.e. closest neighbors first.
            distances : bool
                return the distance of each neighbor to center, along with the
                neighbor itself. Default is False.

        Returns
        -------
            list of Atoms within radius of center or, if distances is True,
            a two-element tuple with (Atom, distance to center).
        """
        if radius <= 0.0:
            raise ValueError(f"Radius must be a positive number: '{radius}'")

        try:
            center_xyz = self.structure.coords[center.index]
        except Exception as err:  # catch all
            if isinstance(center, np.ndarray):
                if center.shape == (3,):
                    center_xyz = center
                else:
                    raise ValueError(
                        f"Center must be a numpy array of shape (1, 3): {center.shape}"
                    )
            else:
                raise TypeError(
                    f"Center must be a valid Atom or numpy array: {type(center)}"
                ) from err

        # No point in turning this into a generator. The C code returns a list.
        nb = self.kdtree.search(center_xyz, radius)

        try:
            idx = center.index
            nb = [p for p in nb if p.index != idx]
        except AttributeError:  # in case of ndarray, no index
            pass

        if sort:
            nb.sort(key=lambda p: p.radius)

        if distances:
            return [(self.structure.topology[p.index], p.radius) for p in nb]
        else:
            return [self.structure.topology[p.index] for p in nb]
