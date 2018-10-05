#!/usr/bin/env python

"""
Analysis of Biomolecular Interfaces.

Module containing interaction type categories and 
analyzers.
"""

from __future__ import print_function

import collections
import logging
import itertools
import os
import warnings

import mdtraj as md 
import numpy as np
import pandas as pd

import simtk.unit as units

from . import functional_groups as fgs

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


# Namedtuple to represent interactions
Interaction = collections.namedtuple('Interaction', ['a', 'b', 'type'])


class InteractionAnalyzerError(Exception):
    """Dummy catch-all class for all exceptions related to `InteractionAnalyzer` objects.
    """
    pass


class InteractionAnalyzer(object):
    """Class to analyze and categorize residue interactions in a Structure object.

    Initialized by providing a parsed structure using one of the OpenMM
    parser classes that provide a `topology` and `positions`.

    Args:
        name (str): path to file used to create Structure instance.
        structure (:obj:`OpenMM Class`): OpenMM `PDB(x)File` object.

    Attributes:
        atom_neighbors (:obj:`list(tuple)`): 2-item tuples with Atom objects.
        neighbors (:obj:`list(tuple)`): 2-item tuples with Residue objects.
    """

    def __init__(self, structure, itable=None):

        self.structure = structure

        if itable is None:
            logging.debug('Creating new InteractionTable')
            _name = 'interactions_{}'.format(os.path.basename(structure.name))
            self.itable = InteractionTable(name=_name)
        else:
            if isinstance(itable, InteractionTable):
                _name = itable.name
                logging.debug('Using existing InteractionTable: {}'.format(_name))
                self.itable = itable
            else:
                emsg = '\'table\' argument should be an InteractionTable() instance'
                raise InteractionAnalyzerError(emsg)

        logging.debug('New InteractionAnalyzer for structure: \'{}\''.format(structure.name))

    #
    # Interaction Calculators
    #

    def _detect_rings(self, residue):
        """Construct spanning tree to determine rings.

        Adapted from https://github.com/ExcitedStates/KGS
            - pdb_structure.py
        """

        def nca_ring(t, v1, v2):
            v1path = []
            v = v1
            while v:
                v1path.append(v)
                v = t[v]
            v = v2
            v2path = []
            while v not in v1path:
                v2path.append(v)
                v = t[v]
            ring = v1path[0:v1path.index(v)+1] + v2path
            return ring

        atoms_in_rings = set()
        parent_dict = {}  # Associates an atom with its parent atom

        neighbors = residue.bonds_per_atom

        for root in residue.atoms():
            if root in parent_dict:
                continue  # Already explored
            parent_dict[root] = None
            fringe = [root]

            while fringe:
                a = fringe[0]
                del fringe[0]

                for n in neighbors[a]:
                    if n in parent_dict and n == parent_dict[a]:
                        continue  # n is just parent of a
                    elif n in parent_dict and not (n in fringe):  # There's a cycle
                        for r in nca_ring(parent_dict, a, n):
                            atoms_in_rings.add(r)
                    elif n not in fringe:
                        parent_dict[n] = a
                        fringe.append(n)

        return atoms_in_rings

    def get_ionic(self, anions=None, cations=None, cutoff=4.0, intramolecular=False, intermolecular=True):
        """Method to detect ionic interactions/salt bridges in the structure.

        Finds potentially charged groups using a list of pre-defined or
        user-supplied functional groups. Functional groups should be defined
        in the `functional_groups` module and be a subclass of `FunctionalGroup`.

        The algorithm expects functional groups divided in anions and cations and
        iteratively searches the structure for these. By default, the groups are:
            cations: `QuaternaryAmine`, `Guanidinium`, `Imidazolium`
            anions: `Carboxylate`, `Phosphate`, `HydrogenPhosphate`, `Sulfate`

        The structure is then searched for pairs of opposing charge groups with
        N/O atoms within, by default, 4 Angstrom of each other (Barlow & Thorton, 1983).

        Args:
            anions (list(:obj: `FunctionalGroup`)): list of anionic groups to match in structure.
            cations (list(:obj: `FunctionalGroup`)): list of cationic groups to match in structure.
            cutoff (float): maximum distance for N/O pairs in matched groups, in Angstrom.
                Default is 4.0 Angstrom
            intramolecular (bool): find ionic interactions *within* the same chain.
                Default is False.
            intermolecular (bool): find ionic interactions only between different chains.
                Default is True.
        """

        # Instantiate ionizable functional groups
        # Anionic
        if anions is None:
            _c = fgs.Carboxylate()
            _p = fgs.Phosphate()
            _ph = fgs.HydrogenPhosphate()
            _s = fgs.Sulfate()

            anions = (_c, _p, _ph, _s)

        # Cationic
        if cations is None:
            _a = fgs.QuaternaryAmine()
            _g = fgs.Guanidinium()
            _i = fgs.Imidazolium()

            cations = (_a, _g, _i)

        # Iterate over residues and find ionizable groups
        res_pos, res_neg = {}, {}
        for res in self.structure.topology.residues():

            for ifg in anions:
                match = ifg.match(res)
                if match:
                    if res not in res_neg:
                        res_neg[res] = set()
                    res_neg[res].update(match)

            for ifg in cations:
                match = ifg.match(res)
                if match:
                    if res not in res_pos:
                        res_pos[res] = set()
                    res_pos[res].update(match)

        logging.debug('Matched {} anionic groups'.format(len(res_neg)))
        logging.debug('Matched {} cationic groups'.format(len(res_pos)))

        # Find pairs that are close in space.
        # Use threshold of 4A between N/O atoms
        _n_or_o = (7, 8)
        _num = 0  # number of interactions for logging
        for res_i, cation in res_pos.items():
            chain_i = res_i.chain.id
            n_and_o = [a for a in cation if a.element.atomic_number in _n_or_o]
            _c, _n, _i, _num_no = chain_i, res_i.name, res_i.id, len(n_and_o)
            logging.debug('Searching neighbors of {}:{}{} ({} N/O)'.format(_c, _n, _i, _num_no))

            neighbors = self.structure.get_neighbors(n_and_o, radius=cutoff, level='atom')

            _seen = set()
            for n in neighbors:
                res_j = n.residue
                chain_j = res_j.chain.id
                if ((res_j not in res_neg) or (res_j in _seen) or
                   (not intramolecular and chain_i == chain_j) or
                   (not intermolecular and chain_i != chain_j)):
                    continue

                j_n_and_o = {a for a in res_neg[res_j] if a.element.atomic_number in _n_or_o}

                if n in j_n_and_o:
                    logging.debug('({})[+] = ({})[-]'.format(res_i, res_j))
                    self.itable.add(res_i, res_j, 'ionic')
                    _num += 1
                    _seen.add(res_j)

        logging.info('Found {} ionic interaction(s) in structure'.format(_num))


class ResidueTable(object):
    """Table-like container object to hold interaction information per-residue.
    """
    pass


class InteractionTable(object):
    """Container class to store and allow search of pairwise interactions.

    Essentially interfaces with a pandas `DataFrame` object, providing utility
    methods to change the content of the `InteractionTable` or query/filter it.
    """

    def __init__(self, name=None):
        self.name = name

        self._setup_new_df()

    def _setup_new_df(self):
        """Creates a new DataFrame.
        """

        _col = ['itype',
                'chain_a', 'chain_b', 'resname_a', 'resname_b',
                'resid_a', 'resid_b']

        self._table = pd.DataFrame(columns=_col)

    def add(self, res_a, res_b, itype=None):
        """Appends an interaction type to the table.

        Args:
            res_a (:obj:`Residue`): `Residue` object involved in the interaction.
            res_b (:obj:`Residue`): Other `Residue` object involved in the interaction.

            itype (str): name of the interaction type (e.g. ionic)
        """

        # Sort residues by chain/number
        res_a, res_b = sorted((res_a, res_b), key=lambda r: (r.chain.id, r.id))

        chain_a, chain_b = res_a.chain.id, res_b.chain.id
        resname_a, resname_b = res_a.name, res_b.name
        resid_a, resid_b = res_a.id, res_b.id

        df = self._table
        df.loc[len(df)] = [itype,
                           chain_a, chain_b, resname_a, resname_b,
                           resid_a, resid_b]

        logging.debug('InteractionTable now contains {} entries'.format(len(df)))

    def clear(self):
        """Deletes all entries in InteractionTable.
        """

        self._table = None
        self._setup_new_df()

    def sort(self):
        pass

    def remove(self):
        pass

    def filter(self):
        pass

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


