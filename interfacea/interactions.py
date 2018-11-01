#!/usr/bin/env python

"""
Analysis of Biomolecular Interfaces.

Module containing interaction type categories and analyzers.
"""

from __future__ import print_function

import collections
import logging
import itertools

import networkx as nx
import numpy as np
import pandas as pd

import simtk.openmm.app.topology as openmm_topology

from . import functional_groups as fgs
from . import data

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())

#
# Classes
#


class InteractionAnalyzerError(Exception):
    """Catch-all class for exceptions related to `InteractionAnalyzer` objects.
    """
    pass


class InteractionAnalyzer(object):
    """Class to analyze and categorize atomic interactions in Structure objects.

    Initialized by providing a `Structure` object.

    Args:
        structure (:obj:`Structure`): interfacea `Structure` object.

    Attributes:
        itable (:obj:`InteractionTable`): queriable table for interactions.
    """

    def __init__(self, structure):

        self.structure = structure

        t_name = 'interactions_{}'.format(structure.name)
        self.itable = InteractionTable(name=t_name)

        # Initialize a few variables that we will need later on
        self.anions = None
        self.aromatics = None
        self.cations = None
        self.hbdonors = None
        self.hydrophobics = None

        msg = 'Created InteractionAnalyzer for structure: \'{}\''
        logging.info(msg.format(structure.name))

    #
    # Private Auxiliary Functions
    #

    def __get_angle(self, a1, a2, a3):
        """Returns the angle between 3 atoms in degrees.

        a1 - a2 - a3 (a2 is the central atom)
        """

        s_xyz = self.structure._np_positions

        vect1 = s_xyz[a2.index] - s_xyz[a1.index]
        vect2 = s_xyz[a2.index] - s_xyz[a3.index]

        v1 = vect1 / np.linalg.norm(vect1)
        v2 = vect2 / np.linalg.norm(vect2)

        # Copied from __vecangle in here to avoid calling self
        return np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

    def __vecangle(self, v1, v2):
        """Returns the angle in degrees between two (normalized) vectors.

        Arguments:
            v1, v2 (:obj:`numpy.array`): N-dimensional numpy arrays of floats.
        """
        return np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

    def __check_planarity(self, atomlist, tolerance=20.0):
        """Returns True if all atoms in atomlist are co-planar.

        Iterates over triplets of atoms (a-b-c) bonded to each other and
        calculates the normal of cross product of the bond vectors (ab, bc).
        If all pairs of normals are parallel (with some margin of error) then
        the atoms are considered planar, else not.

        Arguments:
            atomlist (:obj:`list(Atom)`): list of Atom objects.
            tolerance (float): tolerance factor in degrees to consider a
                deviation from 0/180 still parallel. Defaults to 20.0 degrees.
        """

        tol1, tol2 = 0 + tolerance, 180 - tolerance

        all_xyz = self.structure._np_positions
        vecangle = self.__vecangle

        if len(atomlist) < 3:
            emsg = 'Cannot check planarity for a group of {} atoms: {}'
            raise ValueError(emsg.format(len(atomlist), atomlist))

        normal_list = []
        for atom in atomlist:
            residue = atom.residue
            atom_xyz = all_xyz[atom.index]

            bonded = residue.bonds_per_atom.get(atom)
            for b1, b2 in itertools.product(bonded, repeat=2):

                if b1 == b2:
                    continue

                v1 = all_xyz[b1.index] - atom_xyz
                v2 = all_xyz[b2.index] - atom_xyz

                v1_v2_cross = np.cross(v1, v2)
                v1_v2_cross_u = v1_v2_cross / np.linalg.norm(v1_v2_cross)
                normal_list.append(v1_v2_cross_u)

        for n1, n2 in itertools.product(normal_list, repeat=2):
            angle = vecangle(n1, n2)
            if angle > tol1 and angle < tol2:
                return False

        return True

    def __unpack_subset(self, subset):
        """Unpacks a subset of a structure into a list of residues.

        Arguments:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
        """

        if isinstance(subset, list):
            reslist = []
            for idx, item in enumerate(subset, start=1):
                if isinstance(item, openmm_topology.Residue):
                    reslist.append(item)
                elif isinstance(item, openmm_topology.Chain):
                    for res in item.residues():
                        reslist.append(res)
                else:
                    emsg = 'Item {} (\'{}\') must be a Residue or Chain'
                    raise TypeError(emsg.format(idx, item))
        elif isinstance(subset, openmm_topology.Residue):
            reslist = [subset]
        elif isinstance(subset, openmm_topology.Chain):
            reslist = [res for res in subset.residues()]
        else:
            emsg = '\'{}\' is not a (list of) Residue(s) or Chain(s)'
            raise TypeError(emsg.format(subset))

        return reslist

    #
    # Methods to search and store functional groups
    #

    def find_groups(self, subset=None, group_list=None):
        """Finds functional groups in a structure.

        Arguments:
            group_list (:obj:`list(FunctionalGroup)`): functional groups to
                search the structure for. Must all be instances of the
                `FunctionalGroup` class.
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.

        Returns:
            res_dict (dict): mapping of functional groups to residue.
        """

        if not isinstance(group_list, list):
            emsg = 'Argument \'{}\' must be of type list, not {}'
            raise TypeError(emsg.format(group_list, type(group_list)))

        for idx, item in enumerate(group_list, start=1):
            if not issubclass(item, fgs.FunctionalGroup):
                emsg = 'Item {} (\'{}\') is not a FunctionalGroup subclass'
                raise TypeError(emsg.format(idx, item))

        groups = [fg() for fg in group_list]

        if subset is None:
            reslist = list(self.structure.topology.residues())
        else:
            reslist = self.__unpack_subset(subset)

        # Iterate over residues and find groups
        res_dict = {}
        for res in reslist:
            matches = [set(m) for ifg in groups for m in ifg.match(res)]
            if matches:
                res_dict[res] = matches

        return res_dict

    def find_cations(self, subset=None, cation_list=None):
        """Finds residues matching cationic functional groups.

        Arguments:
            cation_list (:obj:`list(FunctionalGroup)`): positively charged
                groups to search the structure for. Must all be instances of the
                `FunctionalGroup` class. By default, it will search the groups
                defined in `functional_groups.cationic`.
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
        """

        if cation_list is None:
            cation_list = fgs.cationic

        match_dict = self.find_groups(subset=subset, group_list=cation_list)
        n_matches = sum(map(len, match_dict.values()))

        self.cations = match_dict  # cache result
        logging.info('Found {} cationic groups'.format(n_matches))

    def find_anions(self, subset=None, anion_list=None):
        """Finds residues matching anionic functional groups.

        Arguments:
            anion_list (:obj:`list(FunctionalGroup)`): negatively charged
                groups to search the structure for. Must all be instances of the
                `FunctionalGroup` class. By default, it will search the groups
                defined in `functional_groups.anionic`.
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
        """

        # Define & instantiate anions
        if anion_list is None:
            anion_list = fgs.anionic

        match_dict = self.find_groups(subset=subset, group_list=anion_list)
        n_matches = sum(map(len, match_dict.values()))

        self.anions = match_dict  # cache result
        logging.info('Found {} anionic groups'.format(n_matches))

    def find_hydrophobics(self, subset=None, hydrophobic_list=None):
        """Finds residues matching hydrophobic functional groups.

        Arguments:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            hydrophobic_list (:obj:`list(FunctionalGroup)`): hydrophobic groups
                to search the structure for. Must all be instances of the
                `FunctionalGroup` class. By default, it will search the groups
                defined in `functional_groups.anionic`.
        """

        hphobic_list = hydrophobic_list
        if hphobic_list is None:
            hphobic_list = fgs.hydrophobic

        match_dict = self.find_groups(subset=subset, group_list=hphobic_list)

        # Prune sub-matches (e.g. phenyl is in indole)
        for residue, groups in match_dict.items():
            # Sort groups by size, starting with smallest
            sorted_groups = sorted(groups, key=len)
            n_groups = len(sorted_groups)
            for i in range(0, n_groups):
                group_i = sorted_groups[i]
                for j in range(i + 1, n_groups):
                    group_j = set(sorted_groups[j])
                    if group_i.issubset(group_j):
                        msg = 'Residue {}: group {} is a subset of group {}'
                        logging.debug(msg.format(residue, i, j))
                        match_dict[residue].remove(sorted_groups[i])
                        break

        n_matches = sum(map(len, match_dict.values()))

        self.hydrophobics = match_dict  # cache result
        logging.info('Found {} hydrophobic groups'.format(n_matches))

    def find_hb_donors(self, subset=None, donor_list=None):
        """Finds residues matching cationic functional groups.

        Arguments:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            donor_list (:obj:`list(FunctionalGroup)`): list of FunctionalGroup
                object to match hydrogen bond donor groups in the structure. By
                default, it will be the standard `HBondDonor`.
        """

        if donor_list is None:
            donor_list = [fgs.HBondDonor]

        match_dict = self.find_groups(subset=subset, group_list=donor_list)
        n_matches = sum(map(len, match_dict.values()))

        self.hbdonors = match_dict  # cache result
        logging.info('Found {} hydrogen bond donor groups'.format(n_matches))

    def find_aromatic_rings(self, subset=None):
        """Finds residues containing aromatic rings.

        Instead of using functional groups, this method uses a cycle-finding
        function from networkx to find closed cycles in each residue's graph. It
        then checks the planarity of the atoms belonging to the ring and that of
        those directly attached to it.

        Arguments:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.

        Returns:
            a dictionary of AromaticRing objects per residue. An AromaticRing is
            a namedtuple containing the following attributes:
                - residue: residue object containing the ring
                - center: center of mass of the ring
                - radius: maximum distance of any ring atom to center of mass
                - plane: equation defining the ring plane (ax + by + cz + d = 0)
                - normal: one representative of the normal of the ring
        """

        # Define a simple data structure to hold information on each
        # aromatic ring.
        AromaticRing = collections.namedtuple('AromaticRing',
                                              ['residue', 'center', 'radius',
                                               'plane', 'normal'])

        if subset is None:
            reslist = list(self.structure.topology.residues())
        else:
            reslist = self.__unpack_subset(subset)

        res_aromatic = {}
        for residue in reslist:
            atomlist = list(residue.atoms())
            cycle_list = nx.cycle_basis(residue._g)

            # Ensure rings are aromatic by checking planarity
            # For every ring atom, check if the bonded neighbors are in the
            # same plane by comparing normals of the bond vectors
            aromatic_list = []
            for idx, cycle in enumerate(cycle_list):
                cy_atoms = [atomlist[i] for i in cycle]
                is_planar = self.__check_planarity(cy_atoms)
                if is_planar:
                    aromatic_list.append(cy_atoms)

            if not aromatic_list:
                continue

            res_aromatic[residue] = []

            # Calculate properties of the ring: center, radius, plane, normal
            all_xyz = self.structure._np_positions
            for aroring in aromatic_list:
                r_xyz = [all_xyz[a.index] for a in aroring]
                com = np.array(r_xyz).mean(axis=0)  # center of mass
                radius = max((np.linalg.norm(a_xyz - com) for a_xyz in r_xyz))

                # calculate the equation of ring plane
                # Pick 3 points on ring to calculate normal
                p1, p2, p3 = r_xyz[:3]
                v1, v2 = (p2 - p1), (p3 - p1)
                v1_v2_cross = np.cross(v1, v2)
                plane_normal = v1_v2_cross / np.linalg.norm(v1_v2_cross)
                # ax + by + cz = d
                a, b, c = plane_normal
                # Pick atom to get 'd' in equation
                p1x, p1y, p1z = p1
                d = a * p1x + b * p1y + c * p1z
                ring_plane = (a, b, c, d)

                ar = AromaticRing(residue=residue, center=com,
                                  radius=radius, plane=ring_plane,
                                  normal=plane_normal)

                res_aromatic[residue].append(ar)

            n_res_aromatic = len(res_aromatic[residue])
            msg = 'Residue {} contains {} aromatic rings'
            logging.debug(msg.format(residue, n_res_aromatic))

        n_ar_rings = sum(map(len, res_aromatic.values()))
        logging.info('Found {} aromatic rings in structure'.format(n_ar_rings))
        self.aromatics = res_aromatic

    #
    # Interaction Typers/Finders
    #
    def get_ionic(self, subset=None, include_intra=False, max_distance=4.0):
        """Searches structure for interacting groups with opposing charges.

        An interaction between two groups is accepted if any pair of nitrogen or
        oxygen atoms is within a certain distance. By default, the maximum
        distance cutoff is of 4 Angstrom (Barlow & Thorton, 1983).

        Args:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            include_intra (bool): includes intramolecular residue pairs in the
                search. By default, False.
            max_distance (float): maximum distance for N/O pairs in
                matched groups, in Angstrom. Default is 4.0 Angstrom
        """

        max_d = max_distance

        logging.info('Searching structure for interacting charged groups')

        if self.cations is None:
            self.find_cations(subset=subset)
        if self.anions is None:
            self.find_anions(subset=subset)

        cation_dict = self.cations
        anions_dict = self.anions

        # Find pairs of anions/cations close in space
        # Unpack cations nitrogen/oxygen atoms and find neighbors with KDTree
        # Filter results for nitrogen/oxygen of anions.
        def is_nitro_or_oxy(atom):
            return atom.element.atomic_number in {7, 8}

        # Unpack anionic atoms to make it faster to rule out neighbors
        anionic_atoms = {a: (idx_gl, idx_g) for idx_gl, gl in enumerate(anions_dict.values())
                         for idx_g, g in enumerate(gl) for a in g if is_nitro_or_oxy(a)}
        anionic_atoms_set = set(anionic_atoms.keys())

        s = self.structure
        t = self.itable
        _num = 0  # number of interactions for logging
        _seen = set()  # avoid double counting same cat/anion pairs.
        for idx_i, (res, cation_list) in enumerate(cation_dict.items()):
            for idx_ii, cation in enumerate(cation_list):
                cation_no_atoms = filter(is_nitro_or_oxy, cation)

                for at_a in cation_no_atoms:
                    neighbors = s.get_neighbors(at_a, radius=max_d, level='atom')
                    no_list = set(neighbors) & anionic_atoms_set

                    for at_b in no_list:
                        other = at_b.residue

                        if res.chain.id == other.chain.id and not include_intra:
                            continue

                        # idx_j: anion idx, idx_jj: anionic group within anion
                        idx_j, idx_jj = anionic_atoms.get(at_b)
                        pair_id = (idx_i, idx_ii, idx_j, idx_jj)
                        if pair_id in _seen:
                            continue

                        logging.debug('[+] {} - [-] {}'.format(res, other))
                        t.add(res, other, 'ionic', atom_a=at_a, atom_b=at_b)
                        _seen.add(pair_id)
                        _num += 1

        logging.info('Found {} ionic interaction(s) in structure'.format(_num))

    def get_clashes(self, subset=None, include_intra=False, cutoff=1.0):
        """Finds clashes between heavy atoms in the structure.

        Algorithm to find clashes is similar to that used by Chimera's findclash
        tool. Two atoms are considered clashing if their overlap is above a
        certain cutoff value:
            overlap(i, j) = vdw(i) + vdw(j) - d(i,j)
            if overlap(i, j) >= 0.6 the atoms are considered to be clashing.

        This implementation ignores hydrogens, assuming that clashing hydrogens
        can be cleared with a simple minimization and that a clash involving
        heavy atoms is already bad enough. The cutoff is set to 1.0A by default,
        instead of findclash's 0.6, to avoid considering atoms of the peptide
        bond as clashing.

        Args:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            include_intra (bool): includes intramolecular residue pairs in the
                search. By default, False.
            cutoff (float): threshold value in Angstrom to consider that two
                atoms are clashing. It is not a simple distance, see formula
                above. Default is 1.0 Angstrom.
        """

        s = self.structure
        t = self.itable

        def is_hydrogen(atom):
            """Returns True if the atom is a hydrogen.
            """
            return atom.element.atomic_number == 1

        # Determine maximum clash radius
        vdw_radii = data.vdw_radii
        max_radius = max(vdw_radii.values())

        # Search all neighbors within a maximum distance of 2*max_radius
        # of vdw_radius(atom) + max_radius and then use the actual neighbor
        # vdw_radius to screen for clashes.
        clashes = set()
        neighbors = s.get_neighboring_pairs(2.0 * max_radius, level='atom')
        for pair in neighbors:
            atom_i, atom_j, d_ij = pair

            # Ignore same residue and same chain interactions
            res_i, res_j = atom_i.residue, atom_j.residue
            if res_i == res_j:
                continue
            elif (res_i.chain.id == res_j.chain.id) and not include_intra:
                continue

            # Ignore hydrogens
            if is_hydrogen(atom_i) or is_hydrogen(atom_j):
                continue

            # Get vdw(i) + vdw(j)
            vdw_i = vdw_radii.get(atom_i.element.atomic_number)
            vdw_j = vdw_radii.get(atom_j.element.atomic_number)
            overlap = vdw_i + vdw_j - d_ij
            if overlap >= cutoff:

                # Ignore bonded atoms
                if atom_j in res_i.bonds_per_atom[atom_i]:
                    continue

                if (atom_i.residue, atom_j.residue) in clashes:
                    continue

                clashes.add((atom_i.residue, atom_j.residue))
                t.add(res_i, res_j, 'clash', atom_a=atom_i, atom_b=atom_j)

        msg = 'Found {} clashing residues in the structure'
        logging.info(msg.format(len(clashes)))

    def get_hydrophobic(self, subset=None, include_intra=False, max_distance=4.4):
        """Finds interactions between hydrophobic groups.

        A contact is considered between two hydrophobic groups if the minimum
        distance between any pair of non-polar heavy-atoms in the groups is
        below 4.4 Angstrom, as per Table 2 in the reference below:
            Bissantz et al. Journal of Medicinal Chemistry, 2010, vol 53, num 14

        Args:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            include_intra (bool): includes intramolecular residue pairs in the
                search. By default, False.
            max_distance (float): maximum distance for non-polar heavy-atom
                pairs in matched groups, in Angstrom. Default is 4.4 Angstrom
        """

        max_d = max_distance

        logging.info('Searching structure for interacting hydrophobic groups')

        if self.hydrophobics is None:
            self.find_hydrophobics(subset=subset)

        hp_dict = self.hydrophobics

        # Find pairs of hydrophobic groups close in space
        # Unpack apolar heavy-atoms and find neighbors with KDTree
        # Filter results for apolar atoms of hydrophobic groups
        def is_not_h_or_polar(atom):
            return atom.element.atomic_number not in {1, 7, 8}

        _num = 0  # number of interactions for logging
        s = self.structure
        for res, hydrophobic_list in hp_dict.items():

            for idx_i, hp_atoms in enumerate(hydrophobic_list):
                # Find neighbors of atoms in group
                hp_list = list(hp_atoms)
                n_list = s.get_neighbors(hp_list, radius=max_d, level='atom')
                #
                # Make same trick as in cations to save comparisons
                #
                _seen = set()
                for atom in n_list:
                    other = atom.residue
                    if res.chain.id == other.chain.id and not include_intra:
                        continue
                    if other not in hp_dict or not is_not_h_or_polar(atom):
                        continue

                    other_groups = hp_dict[other]
                    for idx_j, group in enumerate(other_groups):
                        pair_id = (idx_i, idx_j)
                        if pair_id in _seen:
                            continue

                        if atom in group:
                            logging.debug('[H] {} - [H] {}'.format(res, other))
                            self.itable.add(res, other, 'hydrophobic')
                            _num += 1
                            _seen.add(pair_id)
                            break

        msg = 'Found {} hydrophobic interaction(s) in structure'
        logging.info(msg.format(_num))

    def get_hbonds(self, subset=None, include_intra=False, max_distance=2.5, min_angle=120.0):
        """Finds hydrogen bonds in the structure using.

        Defines donors as any N/O/F/S connected to a hydrogen atom and
        acceptors as any N/O/F/S within a distance threshold (default 2.5 A) of
        that hydrogen. It then filters the matches for D-H-A triplets with an
        minimum angle threshold (default 120 degrees).

        Similar definitions to those in:
            Baker, E. N., and R. E. Hubbard.
            "Hydrogen bonding in globular proteins."
            Progress in Biophysics and Molecular Biology 44.2 (1984): 97-179

        Args:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            include_intra (bool): includes intramolecular residue pairs in the
                search. By default, False.
            max_distance (float): maximum distance between the hydrogen atom of
                the donor and the acceptor. Default is 4.0 Angstrom
            min_angle (float): minimum angle, in degrees, formed between the
                donor heavy-atom, the donor hydrogen, and the acceptor atom.
                Default is 120 degrees.
        """

        max_d = max_distance
        min_t = min_angle

        # Bring to scope
        get_angle = self.__get_angle
        s = self.structure

        logging.info('Searching structure for hydrogen-bonded interactions')

        # Get donor groups
        if self.hbdonors is None:
            self.find_hb_donors(subset=subset)
        hbdonors_dict = self.hbdonors

        # Find acceptors within range/angle of donors
        _num = 0
        nofs_set = {7, 8, 9, 16}
        occupied = set()  # acceptor indexes
        for res, donor_group in hbdonors_dict.items():

            for group in donor_group:

                # Find acceptors within range
                donor, hydro = group
                acceptors = s.get_neighbors(hydro, radius=max_d, level='atom')

                for acc in acceptors:
                    occupied.add(acc.index)  # acceptor only accepts once
                    other = acc.residue

                    if res.chain.id == other.chain.id and not include_intra:
                        continue
                    elif acc.element.atomic_number not in nofs_set:
                        continue

                    theta = round(get_angle(donor, hydro, acc), 1)  # xxx.y
                    if theta >= min_t:
                        self.itable.add(res, other, 'hbond',
                                        atom_a=donor, atom_b=acc)
                        _num += 1
                        break  # donor can only donate once (X-H pair)

        logging.info('Found {} hydrogen bonds in structure'.format(_num))

    def get_ring_stacking(self, subset=None, include_intra=False, center_d=7.5, angle_dev=30):
        """Searches structure for rings in pi- or t-stacking interactions.

        An interaction between two aromatic rings is accepted if (1) the centers
        of the rings are within a given distance (default 7.5 Angstrom) and (2)
        the planes of the rings form either a 0-30 degree angle (pi-stacking) or
        a 60-90 degree angle (t-stacking). Then, for pi-stacking, the centers of
        the both rings must fall inside the other ring when projected on its
        plane. For t-stacking, both ring centers must be within 5 Angstrom.

        Args:
            subset (:obj:`list(object)`): Residue or Chain objects or list of to
                restrict the search. If None, the default, we will search every
                residue in the structure.
            include_intra (bool): includes intramolecular residue pairs in the
                search. By default, False.
            center_d (float): maximum distance for ring centers, in Angstrom.
            angle_dev (int): maximum deviation, in degrees, for the angles
                between two rings to be considered pi- or t-stacking. Will be
                used to derive a range, e.g.: 0 to angle_dev = pi-stacking and
                90-angle-dev to 90+angle_dev = t-stacking.
        """

        def project_on_plane(point, plane):
            """Finds point on plane nearest to specified point
            """
            x, y, z = point
            a, b, c, d = plane  # normal + d
            t = (d - a * x - b * y - c * z) / (a * a + b * b + c * c)
            return np.array([x + t * a, y + t * b, z + t * c])

        def is_pi_angle(angle, tolerance):
            """Return True if angle meets the pi-stacking criterion
            """
            t1, t2 = 0 + tolerance, 180 - tolerance
            return 0 <= angle <= t1 or t2 <= angle <= 180

        def is_t_angle(angle, tolerance):
            """Return True if angle meets the t-stacking criterion
            """
            return 90 - tolerance <= angle <= 90

        logging.info('Searching structure for ring-stacking interactions')

        if self.aromatics is None:
            self.find_aromatic_rings(subset=subset)

        vecangle = self.__vecangle
        n_pi, n_t = 0, 0
        aromatic_dict = self.aromatics
        reslist = list(aromatic_dict.keys())
        for res_i, res_j in itertools.combinations(reslist, 2):

            if (res_i.chain.id == res_j.chain.id) and not include_intra:
                continue

            rings_i = aromatic_dict[res_i]
            rings_j = aromatic_dict[res_j]
            for ar_i, ar_j in itertools.product(rings_i, rings_j):
                # Distance
                d = np.linalg.norm(ar_i.center - ar_j.center)
                if d > center_d:
                    continue
                # Angle (rounded down to the nearest integer)
                angle = round(vecangle(ar_i.normal, ar_j.normal))
                if is_pi_angle(angle, angle_dev):
                    # Projection Distance
                    proj = project_on_plane(ar_j.center, ar_i.plane)
                    pd = np.linalg.norm(proj - ar_i.center)
                    if pd <= ar_i.radius + 0.5:
                        msg = 'Residues {} & {} are pi-stacking'
                        logging.debug(msg.format(res_i, res_j))
                        self.itable.add(res_i, res_j, itype='pi-stacking')
                        n_pi += 1
                elif is_t_angle(angle, angle_dev):
                    # Stricter Distance
                    if d <= 5.0:
                        msg = 'Residues {} & {} are t-stacking'
                        logging.debug(msg.format(res_i, res_j))
                        self.itable.add(res_i, res_j, itype='t-stacking')
                        n_t += 1

        logging.info('Found {} stacking interactions'.format(n_pi + n_t))


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

        # Identify closest heavy atom in the group and add it as marker
        # of interaction: cation, anion, etc.
        _col = ['itype',
                'chain_a', 'chain_b', 'resname_a', 'resname_b',
                'resid_a', 'resid_b', 'atom_a', 'atom_b', 'energy']

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
        resid_a, resid_b = res_a.id, res_b.id

        atom_a = getattr(kwargs.get('atom_a', None), 'name', None)
        atom_b = getattr(kwargs.get('atom_b', None), 'name', None)

        if _res_a != res_a:  # swapped?
            atom_a, atom_b = atom_b, atom_a

        energy = kwargs.get('energy')

        df = self._table
        df.loc[len(df)] = [itype,
                           chain_a, chain_b, resname_a, resname_b,
                           resid_a, resid_b, atom_a, atom_b, energy]

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
