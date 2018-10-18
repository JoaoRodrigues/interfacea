#!/usr/bin/env python

"""
Analysis of Biomolecular Interfaces.

Module containing interaction type categories and 
analyzers.
"""

from __future__ import print_function

import logging
import os
import warnings

import numpy as np
import pandas as pd

import simtk.unit as units

from . import functional_groups as fgs
from .constants import vdw_radii

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


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
        itable (:obj:`InteractionTable`): queriable table for interactions.
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
    # Auxiliary Functions
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

    def _get_angle(self, a1, a2, a3):
        """Returns the angle between 3 atoms in degrees.

        a1 - a2 - a3 (a2 is the central atom)
        """

        s_xyz = self.structure._np_positions
        vect1 = s_xyz[a2.index] - s_xyz[a1.index]
        vect2 = s_xyz[a2.index] - s_xyz[a3.index]
        l1 = np.sqrt(np.dot(vect1, vect1))
        l2 = np.sqrt(np.dot(vect2, vect2))
        cosa = np.dot(vect1, vect2) / (l1 * l2)
        return np.degrees(np.arccos(cosa))

    #
    # Interaction Typers/Finders
    #

    def get_ionic(self, anions=None, cations=None, cutoff=4.0, intra=False, inter=True):
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
            intra (bool): find ionic interactions *within* the same chain.
                Default is False.
            inter (bool): find ionic interactions only between different chains.
                Default is True.
        """

        logging.info('Finding contacts between charged groups')

        # Instantiate ionizable functional groups
        # Anionic
        if anions is None:
            anions = [fg() for fg in fgs.anionic]
        else:
            if not isinstance(anions, list):
                raise TypeError('Argument \'anions\' must be of type list')
            is_not_fg = [not isinstance(g, fgs.FunctionalGroup) for g in anions]
            if any(is_not_fg):
                emsg = 'All items of \'anions\' must be of type \'FunctionalGroup\''
                raise TypeError(emsg)
            anions = [fg() for fg in anions]

        # Cationic
        if cations is None:
            cations = [fg() for fg in fgs.cationic]
        else:
            if not isinstance(cations, list):
                raise TypeError('Argument \'cations\' must be of type list')
            is_not_fg = [not isinstance(g, fgs.FunctionalGroup) for g in cations]
            if any(is_not_fg):
                emsg = 'All items of \'cations\' must be of type \'FunctionalGroup\''
                raise TypeError(emsg)
            cations = [fg() for fg in cations]

        # Iterate over residues and find ionizable groups
        res_pos, res_neg = {}, {}
        for res in self.structure.topology.residues():
            # returns a list of lists now
            for ifg in anions:
                match_list = ifg.match(res)
                if match_list:
                    if res not in res_neg:
                        res_neg[res] = []
                    res_neg[res].extend(match_list)

            for ifg in cations:
                match_list = ifg.match(res)
                if match_list:
                    if res not in res_pos:
                        res_pos[res] = []
                    res_pos[res].extend(match_list)

        logging.debug('Matched {} anionic groups'.format(len(res_neg)))
        logging.debug('Matched {} cationic groups'.format(len(res_pos)))

        # Find pairs that are close in space.
        # Use threshold of 4A between N/O atoms
        _n_or_o = {7, 8}
        _num = 0  # number of interactions for logging
        for res_i, cation_list in res_pos.items():
            chain_i = res_i.chain.id
            n_and_o = [at for c in cation_list for at in c if at.element.atomic_number in _n_or_o]
            _c, _n, _i, _num_no = chain_i, res_i.name, res_i.id, len(n_and_o)
            logging.debug('Searching neighbors of {}:{}{} ({} N/O)'.format(_c, _n, _i, _num_no))

            neighbors = self.structure.get_neighbors(n_and_o, radius=cutoff, level='atom')

            _seen = set()
            for n in neighbors:
                res_j = n.residue
                chain_j = res_j.chain.id
                if ((res_j not in res_neg) or (res_j in _seen) or
                   (not intra and chain_i == chain_j) or
                   (not inter and chain_i != chain_j)):
                    continue

                anion_list = res_neg[res_j]
                j_n_and_o = {at for a in anion_list for at in a if at.element.atomic_number in _n_or_o}

                if n in j_n_and_o:
                    logging.debug('({})[+] = ({})[-]'.format(res_i, res_j))
                    self.itable.add(res_i, res_j, 'ionic')
                    _num += 1
                    _seen.add(res_j)

        logging.info('Found {} ionic interaction(s) in structure'.format(_num))

    def get_van_der_waals(self):
        pass

    # Add entity keyword to allow narrowing search on a set of atoms
    def get_hydrophobic(self, groups=None, cutoff=4.4, intra=False, inter=True):
        """Finds interactions between hydrophobic groups.

        A contact is considered between two hydrophobic groups if the minimum distance
        between any pair of non-polar heavy-atoms in the groups is below 4.4 Angstrom,
        as per Table 2 in the reference below:
            Bissantz et al. Journal of Medicinal Chemistry, 2010, vol 53, num 14

        Args:
            groups (list(:obj: `FunctionalGroup`)): list of hydrophobic groups to match in
                the structure.
            cutoff (float): maximum distance for non-polar heavy-atom pairs in matched
                groups, in Angstrom. Default is 4.4 Angstrom
            intra (bool): find interactions *within* the same chain.
                Default is False.
            inter (bool): find interactions only between different chains.
                Default is True.
        """

        logging.info('Finding contacts between hydrophobic groups')

        # Define hydrophobic groups
        if groups is None:
            groups = [fg() for fg in fgs.hydrophobic]
        else:
            if not isinstance(groups, list):
                raise TypeError('Argument \'groups\' must be of type list')
            is_not_fg = [not isinstance(g, fgs.FunctionalGroup) for g in groups]
            if any(is_not_fg):
                emsg = 'All items of \'groups\' must be of type \'FunctionalGroup\''
                raise TypeError(emsg)

            groups = [fg() for fg in groups]

        # Match groups to residues in the structure
        matched_residues = {}
        for res in self.structure.topology.residues():
            for ifg in groups:
                match_list = ifg.match(res)
                if match_list:
                    if res not in matched_residues:
                        matched_residues[res] = []
                    matched_residues[res].extend(match_list)

        # Prune matching subsets/supersets (Phenyl & Indole)
        # in the same residue. Keep largest matched group.
        for res in matched_residues:
            to_del = []
            matched_groups = matched_residues[res]
            for i in range(0, len(matched_groups)):
                for j in range(i + 1, len(matched_groups)):
                    if set(matched_groups[i]) <= set(matched_groups[j]):
                        logging.debug('Residue {}: group {} is a subset of group {}'.format(res, i, j))
                        to_del.append(i)
                        break

            matched_residues[res] = [mg for idx, mg in enumerate(matched_groups) if idx not in to_del]

        # Unpack residue to atoms into big set
        logging.debug('Matched {} hydrophobic groups'.format(len(matched_residues)))

        # Calculate distances between non-polar heavy-atoms in groups
        _excl_elem = {1, 7, 8}
        _num = 0  # number of interactions for logging
        for res_i, group_list in matched_residues.items():
            chain_i = res_i.chain.id

            res_npha = [at for c in group_list for at in c if at.element.atomic_number not in _excl_elem]

            _c, _n, _i, _num_npha = chain_i, res_i.name, res_i.id, len(res_npha)
            logging.debug('Searching neighbors of {}:{}{} ({} apolar)'.format(_c, _n, _i, _num_npha))

            neighbors = self.structure.get_neighbors(res_npha, radius=cutoff, level='atom')

            _seen = set()
            for atom_j in neighbors:
                res_j = atom_j.residue
                if res_i == res_j or res_j not in matched_residues or res_j in _seen:
                    continue

                chain_j = res_j.chain.id
                if ((not intra and chain_i == chain_j) or (not inter and chain_i != chain_j)):
                    continue

                logging.debug('({})[+] = ({})[-]'.format(res_i, res_j))
                self.itable.add(res_i, res_j, 'hydrophobic')
                _num += 1
                _seen.add(res_j)

        logging.info('Found {} hydrophobic interaction(s) in structure'.format(_num))

    def get_hydrogen_bonds(self, max_distance=2.5, max_angle=120.0, intra=False, inter=True):
        """Finds hydrogen bonds in the structure.

        Defines acceptors as any N/O/F/S connected to a hydrogen atom and
        donors as any N/O/F/S within 2.5 Angstrom of the hydrogen. It then
        filters the matches for D-H-A triplets with an angle > 120 degrees
        and some heuristics to avoid positively charged nitrogens that do
        not have available lone pairs.

        Similar definitions to those in:
            Baker, E. N., and R. E. Hubbard.
            "Hydrogen bonding in globular proteins."
            Progress in Biophysics and Molecular Biology 44.2 (1984): 97-179
        """

        get_angle = self._get_angle
        s = self.structure

        donor_fg = fgs.HBondDonor()

        _nofs = {7, 8, 9, 16}
        _num = 0
        for res in s.topology.residues():

            donor_chain = res.chain.id
            _seen = set()

            # Find donors
            matches = donor_fg.match(res)
            if not matches:
                continue

            for donor_match in matches:

                # Find acceptors within range
                ha, hydro = donor_match

                _c, _n, _i, _dn = res.chain.id, res.name, res.id, ha.name
                logging.debug('Searching acceptors for {}:{}{}:{}'.format(_c, _n, _i, _dn))
                putative_acceptors = s.get_neighbors(hydro, radius=max_distance, level='atom')

                # Filter acceptors
                # 1. Same residue
                # 2. Chains
                # 3. NOFS
                # 4. Angle
                for a in putative_acceptors:
                    donor_res = a.residue
                    # 1 & 3
                    if donor_res == res or a.element.atomic_number not in _nofs:
                        continue
                    # 2
                    acc_chain = donor_res.chain.id
                    if ((intra and donor_chain != acc_chain) or (inter and donor_chain == acc_chain)):
                        continue
                    # 4
                    theta = get_angle(ha, hydro, a)
                    if theta >= 120.0:
                        self.itable.add(res, donor_res, 'hbond')
                        _num += 1
                        _seen.add(a.residue)  # gotta fix this to allow recording multiple hbonds

        logging.info('Found {} hydrogen bonds in structure'.format(_num))


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


