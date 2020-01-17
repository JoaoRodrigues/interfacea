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

Module containing definitions of functional groups
(subclasses of FunctionalGroup).
"""

from interfacea.chemistry.FunctionalGroup import FunctionalGroup


# Common functional groups

# Charged
class Carboxylate(FunctionalGroup):
    """Carboxylate.

         O
        /
    -- C -- O

    """

    def __init__(self):

        super().__init__(name='carboxylate',
                         charge=-1,
                         elements=[6, 8, 8],
                         bonds=[(0, 1), (0, 2)],
                         max_bonds=[3, 1, 1])


class Carboxyl(FunctionalGroup):
    """Carboxyl.

         O -- H
        /
    -- C -- O

    """

    def __init__(self):
        super().__init__(name='carboxyl',
                         charge=0,
                         elements=[6, 8, 8, 1],
                         # the graph matching algorithm will
                         # match either oxygens to the hydrogen
                         # effectively having 'fuzzy' edges
                         bonds=[(0, 1), (0, 2), (1, 3)],
                         max_bonds=[3, 1, 2, 1])


class Guanidinium(FunctionalGroup):
    """Guanidinium.

          H   H
           \ /
       H    N
       |    |
    -- N -- C
            |
            N
           / \
          H   H
    """

    def __init__(self):
        super().__init__(name='guanidinium',
                         charge=1,
                         elements=[7, 1, 6, 7, 7, 1, 1, 1, 1],
                         bonds=[(0, 1), (0, 2), (2, 3), (2, 4),
                                (3, 5), (3, 6), (4, 7), (4, 8)],
                         max_bonds=[3, 1, 3, 3, 3, 1, 1, 1, 1])


class Imidazole(FunctionalGroup):
    """Imidazole.

    Without protons to avoid ambiguity between ND/NE protonation.
    User can check for position of proton later.
    """

    def __init__(self):
        super().__init__(name='imidazole',
                         charge=0,
                         elements=[6, 6, 7, 6, 7],
                         bonds=[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
                         max_bonds=[3, 3, 3, 3, 3])


class Imidazolium(FunctionalGroup):
    """Imidazolium.
    """

    def __init__(self):
        super().__init__(name='imidazolium',
                         charge=1,
                         elements=[6, 6, 7, 6, 7, 1, 1, 1, 1],
                         bonds=[(0, 1),
                                (1, 2), (1, 5),
                                (2, 3), (2, 6),
                                (3, 4), (3, 7),
                                (4, 0), (4, 8)],
                         max_bonds=[3, 3, 3, 3, 3, 1, 1, 1, 1])


class Phosphate(FunctionalGroup):
    """Phosphate.
    """

    def __init__(self):
        super().__init__(name='phosphate',
                         charge=2,
                         elements=[15, 8, 8, 8],
                         bonds=[(0, 1), (0, 2), (0, 3)],
                         max_bonds=[4, 1, 1, 1])


class SingleCoordinatedPhosphate(FunctionalGroup):
    """Phosphate bound to an atom.
    """

    def __init__(self):
        super().__init__(name='phosphate-h',
                         charge=1,
                         elements=[15, 8, 8, 8, 0],
                         bonds=[(0, 1), (0, 2), (0, 3), (1, 4)],
                         max_bonds=[4, 2, 1, 1, 99])


class QuaternaryAmine(FunctionalGroup):
    """Quaternary Amine.
    """

    def __init__(self):
        super().__init__(name='quaternary_amine',
                         charge=1,
                         elements=[7, 0, 0, 0, 0],
                         bonds=[(0, 1), (0, 2), (0, 3), (0, 4)])


class Sulfonium(FunctionalGroup):
    """Sulfonium.
    """

    def __init__(self):
        super().__init__(name='sulfonium',
                         charge=1,
                         elements=[16, (6, 1), (6, 1), (6, 1)],
                         bonds=[(0, 1), (0, 2), (0, 3)],
                         max_bonds=[3, 4, 4, 4])


class Sulfate(FunctionalGroup):
    """Sulfate.
    """

    def __init__(self):
        super().__init__(name='sulfate',
                         charge=-1,
                         elements=[16, 8, 8, 8],
                         bonds=[(0, 1), (0, 2), (0, 3)],
                         max_bonds=[4, 1, 1, 1])


class HydrogenSulfate(FunctionalGroup):
    """Hydrogen Sulfate.
    """

    def __init__(self):
        super().__init__(name='sulfate-h',
                         charge=0,
                         elements=[16, 8, 8, 8, 1],
                         bonds=[(0, 1), (0, 2), (0, 3), (1, 4)],
                         max_bonds=[4, 2, 1, 1, 1])


# Hydrophobic
class DivalentSulphur(FunctionalGroup):
    """Divalent Sulphur.
    """

    def __init__(self):
        super().__init__(name='divalent-sulphur',
                         charge=0,
                         elements=[16, 0, 0],
                         bonds=[(0, 1), (0, 2)],
                         max_bonds=[2, 4, 4])


class AlkaneCarbon(FunctionalGroup):
    """Alkanes - Aliphatic Saturated Carbon Chain.
    """

    def __init__(self):
        super().__init__(name='alkane',
                         charge=0,
                         elements=[6, (1, 6), (1, 6), (1, 6), (1, 6)],
                         bonds=[(0, 1), (0, 2), (0, 3), (0, 4)])


class AlkeneCarbon(FunctionalGroup):
    """Alkene - Aliphatic Unsaturated Carbon Chain.
    """

    def __init__(self):
        super().__init__(name='alkene',
                         charge=0,
                         elements=[6, (1, 6), (1, 6), (1, 6)],
                         bonds=[(0, 1), (0, 2), (0, 3)],
                         max_bonds=[3, 4, 4, 4])


class AlkyneCarbon(FunctionalGroup):
    """Alkyne - Aliphatic Unsaturated Carbon Chain.
    """

    def __init__(self):
        super().__init__(name='alkene',
                         charge=0,
                         elements=[6, (1, 6), (1, 6)],
                         bonds=[(0, 1), (0, 2), (0, 3)],
                         max_bonds=[2, 4, 4])


class Phenyl(FunctionalGroup):
    """Phenyl Group.
    """

    def __init__(self):
        super().__init__(name='phenyl',
                         charge=0,
                         elements=[6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0],
                         bonds=[(0, 1), (1, 2), (2, 3),
                                (3, 4), (4, 5), (5, 0),
                                (0, 6), (1, 7), (2, 8),
                                (3, 9), (4, 10), (5, 11)])


class Indole(FunctionalGroup):
    """Indole Group.
    """

    def __init__(self):
        super().__init__(name='indole',
                         charge=0,
                         elements=[7, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1, 0, 1, 1, 1, 1],
                         bonds=[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
                                (5, 6), (6, 7), (7, 8), (8, 0), (3, 8),
                                (0, 9), (1, 10), (2, 11), (4, 12), (5, 13),
                                (6, 14), (7, 15)])


class HBondDonor(FunctionalGroup):
    """Hydrogen Bond Donor.
    """

    def __init__(self):
        super().__init__(name='hbond-donor',
                         charge=0,
                         elements=[(7, 8), 1],
                         bonds=[(0, 1)],
                         max_bonds=[99, 1])


# Groups of FunctionalGroups
AnionList = [Carboxylate, Phosphate, SingleCoordinatedPhosphate, Sulfate]
CationList = [Guanidinium, Imidazolium, QuaternaryAmine, Sulfonium]
ChargedList = AnionList + CationList
HydrophobicList = [AlkaneCarbon, Phenyl, Indole, DivalentSulphur]
AromaticList = [Phenyl, Indole]
