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

Module containing interaction type categories and analyzers.
"""

import logging

import pandas as pd

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


class InteractionTable(object):
    """Container class to store and allow search of pairwise interactions.

    Essentially interfaces with a pandas `DataFrame` object, providing utility
    methods to change the content of the `InteractionTable` or query/filter it.
    """

    def __init__(self, name=None):

        self.name = name
        self._setup_new_df()

    def __repr__(self):
        """Pretty print some info.
        """

        df = self._table
        df_str = df.sort_values(by=['itype',
                                    'chain_a', 'chain_b',
                                    'resid_a', 'resid_b']).to_string(index=False)
        return df_str

    def _setup_new_df(self):
        """Creates a new DataFrame.
        """

        # Identify closest heavy atom in the group and add it as marker
        # of interaction: cation, anion, etc.
        _col = ['itype',
                'chain_a', 'chain_b', 'resname_a', 'resname_b',
                'resid_a', 'resid_b', 'atom_a', 'atom_b']

        self._table = pd.DataFrame(columns=_col)

    def add(self, res_a, res_b, itype, **kwargs):
        """Appends an interaction type to the table.

        Arguments:
            res_a (:obj:`Residue`): `Residue` object involved in the interaction.
            res_b (:obj:`Residue`): Other `Residue` object involved in the interaction.
            itype (str): name of the interaction type (e.g. ionic)

        Optional Arguments:
            atom_a (:obj:`Atom`): atom of first residue participating in the
                interaction.
            atom_b (:obj:`Atom`): atom of second residue participating in the
                interaction.
        """

        # Sort residues by chain/number
        _res_a = res_a
        res_a, res_b = sorted((res_a, res_b), key=lambda r: (r.chain.id, r.id))

        chain_a, chain_b = res_a.chain.id, res_b.chain.id
        resname_a, resname_b = res_a.name, res_b.name
        resid_a, resid_b = int(res_a.id), int(res_b.id)

        atom_a = getattr(kwargs.get('atom_a', None), 'name', None)
        atom_b = getattr(kwargs.get('atom_b', None), 'name', None)

        if _res_a != res_a:  # swapped?
            atom_a, atom_b = atom_b, atom_a

        df = self._table
        df.loc[len(df)] = [itype,
                           chain_a, chain_b, resname_a, resname_b,
                           resid_a, resid_b, atom_a, atom_b]

        logging.debug('InteractionTable now contains {} entries'.format(len(df)))

    def clear(self):
        """Deletes all entries in InteractionTable.
        """

        self._table = None
        self._setup_new_df()
        logging.debug('Removed all information from InteractionTable.')

    def to_json(self):
        pass

    def compare(self, other):
        """Performs a per-row/per-column comparison between two tables.

        Proposed algorithm:
            for this.row, other.row in (this, other)
              if this.row not in other
                new.row[itype] = None
              elif this.row in other
                compare this.row[itype] with other.row[itype] => None if diff.
            return new
        """

        pass

    def difference(self, other):
        """Show rows only in this table
        """
        pass
