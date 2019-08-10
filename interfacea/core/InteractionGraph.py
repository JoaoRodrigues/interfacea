#!/usr/bin/env python

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
Analysis of Biomolecular Interfaces.

Module containing the InteractionGraph class for representing molecular
interfaces (and structures) as graphs of residues connected by interactions.
"""

import networkx as nx


class InteractionGraph(object):
    """Data structure to represent molecular interfaces as a queriable graph"""

    def __init__(self, name):
        self.name = name

        self._g = nx.MultiGraph()
        self._node_index = {}  # store residue name
