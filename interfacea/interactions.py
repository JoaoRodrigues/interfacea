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

        self.atom_neighbors = None
        self.neighbors = None

        self.interactions = None  # dictionary of contacts/type

        self.atomic_charges = None
        self.charged_residues = None

        # Convert positions to numpy array
        _xyz_list = structure.positions.value_in_unit(units.angstrom)
        self._positions = np.asarray(_xyz_list, dtype=np.float32)
        del _xyz_list

        logging.debug('New InteractionAnalyzer for structure: \'{}\''.format(structure.name))

    #
    # Neighbor Search Functions
    #
    def get_neighboring_atoms(self, cutoff=5.0):
        """Calculates neighboring atoms based on a minimum distance threshold.

        Args:
            cutoff (float): minimum distance in Angstrom to consider
                two atoms in contact.
        """

        n, m = self._positions.shape
        gram_matrix = np.dot(self._positions, self._positions.T)
        h_matrix = np.tile(np.diag(gram_matrix), (n, 1))
        d_mtx = h_matrix + h_matrix.T - 2 * gram_matrix

        cutoff_sq = cutoff * cutoff
        atom_neighbors = np.argwhere(d_mtx <= cutoff_sq).tolist()

        # Clean list of (i,i) and (i,j)/(j,i) duplicates
        _unique = set()
        for pair in atom_neighbors:
            at_i, at_j = pair
            pair, pair_r = (at_i, at_j), (at_j, at_i)

            if (pair == pair_r) or (pair in _unique) or (pair_r in _unique):
                continue
            _unique.add(pair)

        del atom_neighbors, d_mtx, gram_matrix, h_matrix
        self.atom_neighbors = sorted(_unique)

        logging.debug('Found {} pairs of neighboring atoms at cutoff of {:3.1f}A'.format(len(_unique), cutoff))

    def get_neighboring_residues(self, cutoff=5.0, skip_intra=True):
        """Calculates neighboring residues using minimum atom distances.

        Args:
            cutoff (float): minimum distance in Angstrom to consider
                two atoms in contact and the two residues in contact.

            skip_intra (bool): flag to exclude(default)/include intra-chain
                contacts.

        """

        def _res_sort_key(residue_pair):
            """Private method to sort pairs of residues
            """
            res_i, res_j = residue_pair

            i_key = (res_i.chain.id, res_i.id)
            j_key = (res_j.chain.id, res_j.id)
            return (i_key, j_key)

        if self.atom_neighbors is None:
            self.get_neighboring_atoms(cutoff=cutoff)

        neighbors = {res: [] for res in self.structure.topology.residues()}
        atoms = list(self.structure.topology.atoms())

        for pair in self.atom_neighbors:
            at_i, at_j = pair
            atom_i, atom_j = atoms[at_i], atoms[at_j]

            if skip_intra and (atom_i.residue.chain.id == atom_j.residue.chain.id):
                continue

            if atom_i.residue == atom_j.residue:
                continue

            if atom_j.residue in neighbors[atom_i.residue]:
                continue

            neighbors[atom_i.residue].append(atom_j.residue)
            neighbors[atom_j.residue].append(atom_i.residue)

        self.neighbors = neighbors

        logging.debug('Propagated atom neighbors to residue level')

    #
    # Interaction Calculators
    #

    def _calc_sq_atomic_distance(self, atom_a, atom_b):
        """Returns the squared euclidean distance between two atoms
        """

        pos_diff = self._positions[atom_a.index] - self._positions[atom_b.index]
        return np.sum(pos_diff * pos_diff)

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

    def find_salt_bridges(self, anions=None, cations=None, cutoff=4.0):
        """Method to find salt bridges in the structure.

        Finds potentially charged groups using a list of pre-defined or
        user-supplied functional groups. Functional groups should be defined
        in the `functional_groups` module and be a subclass of `FunctionalGroup`.

        The algorithm expects functional groups divided in anions and cations and
        iteratively searches the structure for these. By default, the groups are:
            cations: `QuaternaryAmine`, `Guanidinium`, `Imidazolium`
            anions: `Carboxylate`, `Phosphate`, `HydrogenPhosphate`, `Sulfate`

        The structure is then searched for pairs of opposing charge groups with
        N/O atoms within, by default, 4 Angstrom of each other (Barlow & Thorton, 1983).
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

        # Run methods if necessary
        if self.neighbors is None:
            self.get_neighboring_residues()

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
        _get_sq_dij = self._calc_sq_atomic_distance
        _sq_cutoff = cutoff * cutoff
        _n_or_o = (7, 8)
        _nsb = 0  # number of salt-bridges
        for res_i, cation in res_pos.items():
            neighbors = self.neighbors[res_i]
            i_n_and_o = [a for a in cation if a.element.atomic_number in _n_or_o]

            for res_j, anion in res_neg.items():
                if res_j not in neighbors:
                    continue

                j_n_and_o = [a for a in anion if a.element.atomic_number in _n_or_o]
                pairs = itertools.product(i_n_and_o, j_n_and_o)
                for ai, aj in pairs:
                    sq_d = _get_sq_dij(ai, aj)
                    logging.debug('({})[+] = ({})[-] == {}'.format(res_i, res_j, sq_d))
                    if sq_d <= _sq_cutoff:
                        self.itable.add(res_i, res_j, 'salt-bridge')
                        _nsb += 1
                        break

        # self.salt_bridges = salt_bridges
        logging.debug('Found {} salt-bridge(s) in structure'.format(_nsb))


class InteractionTable(object):
    """Container class to store and allow search of pairwise interactions.

    Essentially interfaces with a pandas `DataFrame` object, providing utility
    methods to change the content of the `InteractionTable` or query/filter it.
    """

    def __init__(self, name=None):
        self.name = name

        # Setup DataFrame
        _col = ['itype',
                'chain_a', 'chain_b', 'resname_a', 'resname_b',
                'resid_a', 'resid_b', 'energy']

        self._table = pd.DataFrame(columns=_col)

    def add(self, res_a, res_b, itype=None, ienergy=None):
        """Appends an interaction type to the table.

        Args:
            res_a (:obj:`Residue`): `Residue` object involved in the interaction.
            res_b (:obj:`Residue`): Other `Residue` object involved in the interaction.

            itype (str): name of the interaction type (e.g. salt-bridge)
            ienergy (float): energy value of the interaction.
        """

        # Sort residues by chain/number
        res_a, res_b = sorted((res_a, res_b), key=lambda r: (r.chain.id, r.id))

        chain_a, chain_b = res_a.chain.id, res_b.chain.id
        resname_a, resname_b = res_a.name, res_b.name
        resid_a, resid_b = res_a.id, res_b.id

        df = self._table
        df.loc[len(df)] = [itype,
                           chain_a, chain_b, resname_a, resname_b,
                           resid_a, resid_b, ienergy]

        logging.debug('InteractionTable now contains {} entries'.format(len(df)))

    def sort(self):
        pass

    def remove(self):
        pass

    def filter(self):
        pass

    def to_json(self):
        pass