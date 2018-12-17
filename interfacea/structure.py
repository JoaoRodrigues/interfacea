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

Module containing Structure class to represent 3D molecular objects and allow
fundamental manipulation/parameterization functions.
"""

import itertools
import logging
import os
import tempfile
import sys
import warnings

import numpy as np
import networkx as nx

import pdbfixer as pf
from pdbfixer.pdbfixer import Sequence

import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as units

from . import data
from .src import kdtrees
from .private import internal

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


class StructureError(Exception):
    """Dummy catch-all class for all exceptions related to `Structure` objects.
    """
    pass


class Structure(object):
    """Wrapper around OpenMM classes to store/manipulate structural data.

    Initialized by providing a parsed structure using one of the OpenMM
    parser classes that provide a `topology` and `positions`.

    Args:
        name (str): path to file used to create Structure instance.
        structure (:obj:`OpenMM Class`): OpenMM `PDB(x)File` object.
        build_kdtree (bool): automatically build KDTree on instantiation.
            Default is True.

    Attributes:
        topology (:obj:`OpenMM Topology`): OpenMM topology.
        positions (:obj:`OpenMM Positions`): OpenMM positions array.

        sequences (:obj:`list(`PDBFixer Sequence`)`: list of sequences described
            in the `Structure` object.
        forcefield (str): name of the forcefield used to parameterize the
            structure. Default is None - the structure is not parameterized.
        potential_energy (float): potential energy (in kJ/mol) calculated using
            the forcefield parameters.
    """

    def __init__(self, name, structure, build_kdtree=True):

        self.sequences = None
        self.forcefield = None  # forcefield name (str)
        self.potential_energy = None

        self._forcefield = None  # forcefield object
        self._system = None
        self._pdbfixer = None  # cache PDBFixer structure if we need it.
        self._kdt = None

        self.name = name
        self.__set_positions(structure.positions)
        self.__set_topology(structure.topology)

        # Build KDTree
        if build_kdtree:
            self.__build_kdtree()

        logging.info('Created Structure from \'{}\''.format(name))

    def __repr__(self):
        """Print pretty things when called.
        """

        rep_str = 'Structure ({})'.format(os.path.basename(self.name))

        n_chain = self.topology.getNumChains()
        n_resid = self.topology.getNumResidues()
        n_atoms = self.topology.getNumAtoms()

        rep_str += ' ({} chain(s), '.format(n_chain)
        rep_str += '{} residue(s), '.format(n_resid)
        rep_str += 'and {} atom(s))'.format(n_atoms)

        if self._forcefield:
            rep_str += ' (ff={})'.format(self.forcefield)

        return rep_str

    # Pickling/Copying methods
    def copy(self):
        """Returns a (deep) copy of the Structure.

        Does not include KDTree or System objects to allow pickling.
        """
        newstruct = self.__class__(self.name, self, build_kdtree=False)
        warnings.warn('Copy is not parameterized.')  # kdtree is built on search
        return newstruct

    def __copy__(self, *args):
        """Shallow Copy Override
        """

        c = self.copy()
        return c

    def __deepcopy__(self, *args):
        """Deepcopy override
        """

        c = self.copy()
        return c

    #
    # Private Methods
    #

    def __build_kdtree(self):
        """Builds a KDTree for fast neighbor lookup.

        Biopython's KDTree in C and Python bindings by Michiel de Hoon.
        For details, read the source code at src/kdtrees.c

        The module implements a KDTree class that takes a Nx3 numpy array of
        doubles.

        The class exposes three methods:
            search: returns indices of atoms within a radius (in Angstrom)
                    around a central point.
            neighbor_search: returns all point pairs within a certain radius of
                             each other.
            neighbor_search_simple: same as above, slow implementation for test
                                    purposes.

        In interfacea, we will expose search and neighbor_search but wrap them
        to return either Atom, Residue, or Chain objects, a bit like Biopython.
        """

        logging.debug('Building KDTree (this might take a minute or two)')
        self._kdt = kdtrees.KDTree(self._np_positions)

    def __set_topology(self, topology):
        """Utility method to apply changes on topology changes.
        """

        logging.debug('(Re)Setting topology for \'{}\''.format(self.name))
        self.topology = topology

        self.__get_bonded_atoms()
        self.__make_residue_graphs()

    def __set_positions(self, positions):
        """Utility method to apply changes on atom addition/deletion.
        """

        logging.debug('(Re)Setting coordinates for \'{}\''.format(self.name))
        self.positions = positions

        # Convert positions to numpy array
        # use double precision for coordinates because of KDTree
        _xyz_list = positions.value_in_unit(units.angstrom)
        self._np_positions = np.asarray(_xyz_list, dtype="d")
        del _xyz_list

        # Update KDTree
        if self._kdt is not None:
            self.__build_kdtree()

    def __load_to_pdbfixer(self):
        """Class to write a temporary PDB file and reload using PDBFixer.

        If PDBFixer was never called before, runs and caches the resulting
        Structure object. Always resets/empties the missing lists to avoid
        conflicts.
        """

        if self._pdbfixer is None:
            with tempfile.TemporaryFile(mode='r+') as handle:
                app.PDBFile.writeFile(self.topology, self.positions, handle,
                                      keepIds=True)
                handle.seek(0)  # rewind
                s = pf.PDBFixer(pdbfile=handle)

            self._pdbfixer = s

        sequences = []
        for chain in self._pdbfixer.topology.chains():
            chain_reslist = [r.name for r in chain.residues()]
            sequences.append(Sequence(chain.id, chain_reslist))

        self._pdbfixer.sequences = self.sequences = sequences
        self._pdbfixer.missingAtoms = {}
        self._pdbfixer.missingResidues = {}
        self._pdbfixer.missingTerminals = {}

        logging.debug('Cached PDBFixer Data Structure')

    def __guess_bonds_from_coordinates(self, residue):
        """Guesses atom connectivity from atomic distances and covalent radii.

        This is only meant to be used in the absence of proper parameters in the
        forcefield for a particular residue. Will patch, hopefully, enough info
        for the graph matching to work, but forget about minimizations/energies.

        Args:
            residue (:obj:`Residue`): residue object to define connectvity for.
                Populates the `Residue.bonds()` attribute in-place.
        """

        # Taken from m4xbondage.html webpage
        metals = {
            3, 4, 12, 13, 18, 23, 24, 25, 27, 29, 30, 31, 33, 36, 37, 38, 39,
            42, 47, 48, 49, 51, 52, 54, 55, 56, 57, 58, 62, 63, 64, 65, 67, 70,
            71, 74, 75, 76, 77, 78, 79, 80, 81, 82, 92,
        }

        tolerance = 0.45  # tolerance value to account for weird structures

        xyz = self._np_positions
        cov_radii = data.covalent_radii

        msg = 'Determining atom connectivity from coordinates for residue: {}'
        logging.info(msg.format(residue.name))

        _num = 0  # number of bonds found
        atomlist = list(residue.atoms())
        for atom_i, atom_j in itertools.combinations(atomlist, 2):

            elem_i = atom_i.element.atomic_number
            elem_j = atom_j.element.atomic_number
            if elem_i in metals or elem_j in metals:
                continue  # skip bonds with metallic elements

            radius_i = cov_radii.get(elem_i, 2.0)
            radius_j = cov_radii.get(elem_j, 2.0)

            d_ij = np.linalg.norm(xyz[atom_i.index] - xyz[atom_j.index])
            if d_ij <= (radius_i + radius_j + tolerance):
                self.topology.addBond(atom_i, atom_j)
                _num += 1

        msg = 'Assigned {} bonds in residue {}'
        logging.debug(msg.format(_num, residue.name))

    def __get_bonded_atoms(self):
        """Creates a dictionary of bonds per atom.

        Avoids residue.bonds() purposedly because it is very slow.
        """

        logging.debug('__get_bonded_atoms :: started')

        # Make a nested dictionary of residue bonds {res: atom: [bonds]}
        top = self.topology
        res_bonds = {r: {a: [] for a in r.atoms()} for r in top.residues()}

        # Iterate over topology bonds
        for a1, a2 in self.topology.bonds():
            r1, r2 = a1.residue, a2.residue
            # Skip bonds between residues
            if r1 != r2:
                continue

            res_bonds[r1][a1].append(a2)
            res_bonds[r2][a2].append(a1)

        # Iterate over each residue and assign to attribute
        for res in self.topology.residues():
            bond_dict = res_bonds[res]
            n_bonds = sum(map(len, bond_dict.values()))

            if not n_bonds:  # unknown residue?
                self.__guess_bonds_from_coordinates(res)
                # Do not except many, so OK to use .internal_bonds() here
                for a1, a2 in res.internal_bonds():
                    res_bonds[res][a1].append(a2)
                    res_bonds[res][a2].append(a1)

                bond_dict = res_bonds[res]
                n_bonds = sum(map(len, bond_dict.values()))
                if not n_bonds:
                    wmsg = 'Residue {}:{}{} is missing bonding information.'
                    warnings.warn(wmsg.format(res.chain.id, res.name, res.id))

            res.bonds_per_atom = bond_dict

        logging.debug('__get_bonded_atoms :: finished')

    def __make_residue_graphs(self):
        """Creates a networkx.Graph representation of a `Residue`.

        Uses the atom elements as node attributes and the topology bonds
        as edges.

        Somewhat a performance bottleneck for large structures.
        Should we optimize? How?
        """

        logging.debug('__make_residue_graphs :: started')

        reslist = list(self.topology.residues())
        for residue in reslist:
            at_to_idx = {at: idx for idx, at in enumerate(residue.atoms())}

            # Make graph of residue
            res_g = nx.Graph()
            for atom, idx in at_to_idx.items():
                res_g.add_node(idx, element=atom.element.atomic_number)
                for bonded in residue.bonds_per_atom.get(atom, []):
                    res_g.add_edge(idx, at_to_idx[bonded])

            residue._g = res_g

        logging.debug('__make_residue_graphs :: finished')

    def __load_forcefield(self, forcefield='amber14-all.xml'):
        """Utility private method to load forcefield definitions.
        """

        try:
            loaded_forcefield = app.ForceField(forcefield)
            logging.debug('Loaded forcefield: {}'.format(forcefield))

        except ValueError as e:
            emsg = 'Error when loading forcefield XML file: {}'
            raise StructureError(emsg.format(forcefield)) from e

        self._forcefield = loaded_forcefield
        self.forcefield = forcefield

    #
    # Public Methods
    #

    # IO
    def write(self, output, ftype=None, overwrite=False):
        """Writes `Structure` object to file.

        Uses OpenMM PDBFile or PDBxFile methods to write the `Structure` to a
        file on disk in PDB or mmCIF format, respectively. The output format is
        guessed from the user-provided file name or by the optional argument
        `format`.


        Args:
            output (str): name of file to write Structure to disk.
            ftype (str): file format to use when writing the file. Must be
                either 'pdb' or 'cif'. If None, tries guessing from output file
                extension.
            overwrite(bool, optional): write file even if it already exists. Defaults to False.

        Raises:
            StructureError: if file type is not supported.
            OSError: if file already exists and overwrite is set to False.
        """

        _writers = {'cif': app.PDBxFile.writeFile, 'pdb': app.PDBFile.writeFile}
        _fmt_str = ','.join(_writers.keys())

        if ftype is None:  # read from filename
            _, ext = os.path.splitext(output)
            ftype = ext[1:]  # removes the dot
            if not ftype.strip():  # empty (no extension?)
                emsg = 'File type could not be guessed from output name: {}'
                raise StructureError(emsg.format(output))

        writer = _writers.get(ftype)
        if writer is None:
            emsg = 'Unsupported file type \'{}\'. Choose from {}'
            raise StructureError(emsg.format(ftype, _fmt_str))

        if isinstance(output, str):
            if os.path.isfile(output) and not overwrite:
                emsg = 'File already exists. Use overwrite=True or remove file.'
                raise OSError(emsg)
            try:
                with open(output, 'w') as handle:
                    writer(self.topology, self.positions, handle, keepIds=True)
            except Exception as e:
                emsg = 'Error when writing Structure to file: {}'
                raise StructureError(emsg.format(handle.name)) from e
        else:
            raise TypeError('\'output\' argument must be a string.')

    # Structure manipulation
    def remove_solvent(self):
        """Removes solvent molecules from the structure.
        """

        solvent = set(('HOH', 'WAT'))

        m = app.Modeller(self.topology, self.positions)
        solvent_list = [r for r in m.topology.residues() if r.name in solvent]
        m.delete(solvent_list)

        msg = 'Removed {} solvent molecules from Structure'
        logging.info(msg.format(len(solvent_list)))

        self.__set_topology(m.topology)
        self.__set_positions(m.positions)

    def prepare(self, cap_termini=True, forcefield='amber14-all.xml', minimize=True, pH=7.0):
        """Utility function to complete a structure and minimize it.

        Essentially performs the same as calling the following:
            mol = interfacea.read('...')
            mol.add_termini()
            mol.add_missing_atoms()
            mol.protonate(pH=7.0)
            mol.parameterize()
            mol.minimize(iterations=250, hydrogen_only=False)
        """

        if cap_termini:
            self.add_termini()

        self.add_missing_atoms()
        self.protonate(pH=pH)
        self.parameterize(forcefield=forcefield)

        if minimize:
            self.minimize()

    def add_termini(self, ends=None):
        """Adds missing terminal atoms to protein chains in the `Structure`.

        Uses PDBFixer to modify structure. By default, adds acetyl (ACE) and
        N-methyl (NME) caps to N- and C- termini of protein molecules. Other
        caps can be specified using the `ends` option, as long as they are
        supported by PDBFixer and OpenMM.

        Args:
            ends (:obj:`list(tuple(str, str)`, optional): definitions of termini
                groups to add to each chain. Number of items in list must match
                the number of chains in the `Structure`. Allowed options for
                each chain are: 'ACE', 'NME', or None (charged terminus).

        Raises:
            StructureError
        """

        protein_aa = data.protein_aa

        _allowed_n_caps = set(('ACE', None))
        _allowed_c_caps = set(('NME', None))

        chains = list(self.topology.chains())
        num_chains = len(chains)
        logging.debug('Adding termini caps to {} chains'.format(num_chains))

        if ends is not None:
            num_ends = len(ends)
            if num_ends != num_chains:
                emsg = 'Number of terminal capping groups ({}) does not match '
                emsg += 'the number of chains in the molecule ({})'
                raise StructureError(emsg.format(num_ends, num_chains))

            for chain_idx, chain in enumerate(ends):
                name = chains[chain_idx].id
                n_cap, c_cap = chain

                if n_cap not in _allowed_n_caps:
                    emsg = 'Unsupported N-terminal cap for chain {}: {}'
                    raise StructureError(emsg.format(name, n_cap))

                if c_cap not in _allowed_c_caps:
                    emsg = 'Unsupported C-terminal cap for chain {}: {}'
                    raise StructureError(emsg.format(name, c_cap))

        else:
            ends = [('ACE', 'NME') for _ in range(num_chains)]

        # PDBFixer is picky with mmCIF files and OpenMM does not
        # export them properly. The irony I know.
        # For now, we save as PDB and re-read.
        self.__load_to_pdbfixer()
        s = self._pdbfixer

        # Find missing atoms to exclude from termini addition
        # Avoids adding other atoms when calling `add_termini`
        # Clunky but explicit is better than implicit!
        s.findMissingResidues()
        s.findMissingAtoms()
        missing_residues = set(s.missingResidues.keys())
        missing_atoms = set(s.missingAtoms.keys())

        # Build Sequence object(s) including caps
        sequences = []
        for chain_idx, chain in enumerate(s.topology.chains()):
            chain_reslist = [r.name for r in chain.residues()]
            n_ter, c_ter = chain_reslist[0], chain_reslist[-1]

            # Add caps if necessary and if protein molecules
            n_cap, c_cap = ends[chain_idx]
            if n_cap and n_ter != n_cap and n_ter in protein_aa:
                chain_reslist.insert(0, n_cap)
                msg = 'Adding \'{}\' capping group to chain {} N-terminus'
                logging.debug(msg.format(n_cap, chain.id))
            else:
                msg = 'Ignoring N-terminal caps on chain {}'
                logging.debug(msg.format(chain.id))

            if c_cap and c_ter != c_cap and c_ter in protein_aa:
                chain_reslist.append(c_cap)
                msg = 'Adding \'{}\' capping group to chain {} C-terminus'
                logging.debug(msg.format(c_cap, chain.id))
            else:
                msg = 'Ignoring C-terminal caps on chain {}'
                logging.debug(msg.format(chain.id))

            sequences.append(Sequence(chain.id, chain_reslist))

        s.sequences = sequences
        self._pdbfixer.sequences = self.sequences = sequences

        # Resume filtering of missing atoms pertaining to non-terminal residues.
        s.findMissingResidues()
        for res in list(s.missingResidues.keys()):
            if res in missing_residues:
                del s.missingResidues[res]

        s.findMissingAtoms()
        for res in list(s.missingAtoms.keys()):
            if res in missing_atoms:
                del s.missingAtoms[res]

        # Add missing terminal atoms/residues
        try:
            s.addMissingAtoms(seed=internal.RANDOM_SEED)
        except TypeError as err:
            msg = 'OpenMM does not add missing atoms at reproducible positions.'
            msg += ' Consider installing the latest OpenMM from git'
            msg += ' to minimize this effect. This is just a warning, nothing'
            msg += ' is wrong with your structure or the added atoms.'
            warnings.warn(msg)

            s.addMissingAtoms()

        self.__set_topology(s.topology)
        self.__set_positions(s.positions)

    def add_missing_atoms(self):
        """Method to add missing atoms to a `Structure` object.

        Uses PDBFixer to analyze structure and add missing heavy atoms.
        """

        # Save as PDB and re-read.
        # Someone should really write an OpenMM to PDBFixer conversion ...
        self.__load_to_pdbfixer()
        s = self._pdbfixer

        s.findMissingResidues()
        s.findMissingAtoms()

        s.missingTerminals = {}  # do not add missing terminals here.

        n_added_atoms = len(s.missingAtoms)

        if n_added_atoms:
            logging.info('Found missing heavy atoms: {}'.format(n_added_atoms))

            try:
                s.addMissingAtoms(seed=internal.RANDOM_SEED)
            except TypeError as err:
                msg = 'OpenMM does not add missing atoms at reproducible positions.'
                msg += ' Consider installing the latest OpenMM from git'
                msg += ' to minimize this effect. This is just a warning, nothing'
                msg += ' is wrong with your structure or the added atoms.'
                warnings.warn(msg)

                s.addMissingAtoms()

            self.__set_topology(s.topology)
            self.__set_positions(s.positions)

            # Issue warning about atom positions
            wmsg = ('Atoms added but their positions are not optimized. '
                    'Minimize the structure before doing any analysis')
            warnings.warn(wmsg)

    def protonate(self, forcefield='amber14-all.xml', pH=7.0, keep_existing=False):
        """Method to add hydrogen atoms to a `Structure` object.

        Uses the `Modeller` class from OpenMM to add hydrogen atoms to the
        structure. Removes existing hydrogen atoms (to avoid naming issues)
        before adding them again with naming and topology matching the
        forcefield and chosen pH.


        Args:
            forcefield (str): name of file defining the forcefield.
            pH (float): numerical value of the pH to check protonation states
                of ionizable groups.
            keep_existing (bool): Does not remove existing protons. This can be
                dangerous because protons might not follow proper naming
                conventions but useful in the case of wanting specific
                protonation states. Default is False.
        """

        def is_hydrogen(atom):
            """Returns True if the atom is a hydrogen atom.
            """
            return atom.element.atomic_number == 1

        if self._forcefield is None or self.forcefield != forcefield:
            self.__load_forcefield(forcefield)

        model = app.Modeller(self.topology, self.positions)

        if not keep_existing:
            existing_H = [a for a in model.topology.atoms() if is_hydrogen(a)]
            model.delete(existing_H)
            msg = 'Removed {} existing hydrogen atoms'
            logging.debug(msg.format(len(existing_H)))

        logging.debug('Protonating structure at pH {}'.format(pH))
        model.addHydrogens(forcefield=self._forcefield, pH=pH)

        self.__set_topology(model.topology)
        self.__set_positions(model.positions)

        # Issue warning about atom positions
        wmsg = ('Protons added but their positions are not optimized. '
                'Minimize the structure before doing any analysis')
        warnings.warn(wmsg)

    def mutate(self, mutation_list):
        """Mutates residues in the molecule using PDBFixer.

        This is a very crude method of deleting/adding atoms, so most useful
        (or reasonable) for single mutations. Mutation of multiple residues at
        once might yield a very bad structure, even if followed by minimization.

        Args:
            mutation_list (:obj:`list(tuple)`): list of two-item tuples
                containing the id of the residue to mutate as a string with
                'chain-resname-resid' and the resname of the mutated residue.
                Residue names should always be in three-letter code to avoid
                ambiguities:
                    e.g. ('A-ASN-1', 'ALA') mutates ASN1 of chain A to alanine.

        Raises:
            StructureError
        """

        # Build list of valid residues to mutate to
        _pf = pf.pdbfixer
        resnames = _pf.proteinResidues + _pf.dnaResidues + _pf.rnaResidues
        _supported_resnames = set(resnames)

        # Sanity check on mutation list
        mut_per_chain = {}

        if not isinstance(mutation_list, list):
            if isinstance(mutation_list, tuple):
                mutation_list = [mutation_list]  # assume single mutation
            else:
                emsg = '\'mutation_list\' must be a list. Check documentation.'
                raise TypeError(emsg)

        for mutation in mutation_list:
            try:
                ori, mutres = mutation
                chain, resname, idx = ori.split('-')
            except Exception as e:
                emsg = 'Wrong format in mutation: \'{}\''.format(mutation)
                raise StructureError(emsg) from e

            if mutres not in _supported_resnames:
                emsg = 'Residue not supported for mutation: {}'.format(mutres)
                raise StructureError(emsg)

            # defer to PDBFixer to catch errors of mutating non-existing
            # residues or on non-existing chains.
            if chain not in mut_per_chain:
                mut_per_chain[chain] = []

            msg = 'Mutating residue {}:{}{} to {}'
            logging.info(msg.format(chain, resname, idx, mutres))
            mut_per_chain[chain].append('{}-{}-{}'.format(resname, idx, mutres))

        self.__load_to_pdbfixer()
        s = self._pdbfixer

        # Mutate only. Do not add/complete structure.
        s.findMissingResidues()
        s.findMissingAtoms()
        missing_residues = set(s.missingResidues.keys())
        missing_atoms = set(s.missingAtoms.keys())

        # Mutate on each chain at a time
        for chain in mut_per_chain:
            muts = mut_per_chain[chain]

            try:
                s.applyMutations(muts, chain)
            except (KeyError, ValueError) as e:
                emsg = 'Unknown error when applying mutations to the structure'
                raise StructureError(emsg) from e

            s.findMissingResidues()
            for res in list(s.missingResidues.keys()):
                if res in missing_residues:
                    del s.missingResidues[res]

            s.findMissingAtoms()
            for res in list(s.missingAtoms.keys()):
                if res in missing_atoms:
                    del s.missingAtoms[res]

            s.missingTerminals = {}

            try:
                s.addMissingAtoms(seed=internal.RANDOM_SEED)
            except TypeError as err:
                msg = 'OpenMM does not add missing atoms at reproducible positions.'
                msg += ' Consider installing the latest OpenMM from git'
                msg += ' to minimize this effect. This is just a warning, nothing'
                msg += ' is wrong with your structure or the added atoms.'
                warnings.warn(msg)

                s.addMissingAtoms()

        self.__set_topology(s.topology)
        self.__set_positions(s.positions)

        # Issue warning about atom positions
        wmsg = ('Residue mutated but atom positions were not optimized. '
                'Minimize the structure before doing any analysis')
        warnings.warn(wmsg)

    #
    # MM Functions
    #
    def parameterize(self, forcefield='amber14-all.xml'):
        """Wrapper function to create an OpenMM system from the Structure.

        By default, does not use a cutoff (nor I plan to change this).

        Args:
            forcefield (str): name of 'xml' file containing forcefield definitions.
        """

        if self._forcefield is None or forcefield != self.forcefield:
            if forcefield != self.forcefield:
                wmsg = 'Structure previously parameterized with \'{}\'. '
                wmsg += 'It is advisable to run protonate() again'
                warnings.warn(wmsg.format(self.forcefield))

            self.__load_forcefield(forcefield)

        ff = self._forcefield
        system = ff.createSystem(self.topology, nonbondedMethod=app.NoCutoff)
        self._system = system

        logging.debug('Structure parameterized using \'{}\''.format(forcefield))

    def calculate_energy(self):
        """Calculates the potential energy of the system.
        """

        if self._system is None:
            emsg = 'Structure is not parameterized yet. Use `parameterize()`.'
            raise StructureError(emsg)

        # Set integrator
        integrator = mm.LangevinIntegrator(300 * units.kelvin,
                                           1.0 / units.picosecond,
                                           2.0 * units.femtosecond)
        integrator.setRandomNumberSeed(internal.RANDOM_SEED)
        integrator.setConstraintTolerance(0.00001)

        # Create context
        context = mm.Context(self._system, integrator)
        context.setPositions(self.positions)

        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        energy_kjmol = energy.value_in_unit(units.kilojoule_per_mole)
        msg = 'Potential energy of the structure = {:8.3f} kJ/mol'
        logging.info(msg.format(energy_kjmol))

        if energy_kjmol > 0.0:  # sort of arbitrary value
            wmsg = 'Potential energy of the structure is high ({:8.3f}). '
            wmsg += 'You should try minimize() before any analysis'
            warnings.warn(wmsg.format(energy_kjmol))

        self.potential_energy = energy

        del context
        del integrator

    def minimize(self, iterations=250, hydrogen_only=False):
        """Perform energy minimization using OpenMM.

        Args:
            iterations (int): number of steps taken by the minimizer
                (or 0 until convergence). Default is 250 steps.
            hydrogen_only (bool): minimize positions of hydrogen atoms only.
                Default is False - all atoms are free to move.
        """

        if self._system is None:
            emsg = 'Structure is not parameterized yet. Use `parameterize()`.'
            raise StructureError(emsg)

        # Set integrator
        integrator = mm.LangevinIntegrator(300 * units.kelvin,
                                           1.0 / units.picosecond,
                                           2.0 * units.femtosecond)
        integrator.setRandomNumberSeed(internal.RANDOM_SEED)
        integrator.setConstraintTolerance(0.00001)

        # Harmonic position restraints
        posre_k = 500.0 * (units.kilojoule_per_mole / units.nanometer**2)
        posre = mm.CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
        posre.addGlobalParameter("k", posre_k)
        posre.addPerParticleParameter("x0")
        posre.addPerParticleParameter("y0")
        posre.addPerParticleParameter("z0")

        if hydrogen_only:

            logging.info('Adding position restraints to all heavy atoms')
            elemlist = [a.element.atomic_number for a in self.topology.atoms()]

            for idx, xyz in enumerate(self.positions):
                elem = elemlist[idx]
                if elem != 1:
                    posre.addParticle(idx, xyz.value_in_unit(units.nanometers))

            self._system.addForce(posre)

            n_restraints = posre.getNumParticles()
            n_atoms = len(elemlist)
            msg = 'Restrained {} out of {} atoms'
            logging.debug(msg.format(n_restraints, n_atoms))

        # Perform minimization
        context = mm.Context(self._system, integrator)
        context.setPositions(self.positions)

        state = context.getState(getEnergy=True)
        initial_e = state.getPotentialEnergy()
        initial_e_kjmol = initial_e.value_in_unit(units.kilojoule_per_mole)
        msg = 'Energy before minimization: {:8.3f} kJ/mol'
        logging.info(msg.format(initial_e_kjmol))

        if initial_e_kjmol > 1000000.0:  # sort of arbitrary value
            wmsg = ('Initial potential energy is very high ({:8.3f} kJ/mol). '
                    'Minimization is likely to fail.'.format(initial_e_kjmol))
            warnings.warn(wmsg)

        mm.LocalEnergyMinimizer.minimize(context, maxIterations=iterations)

        state = context.getState(getPositions=True, getEnergy=True)
        final_e = state.getPotentialEnergy()
        final_e_kjmol = final_e.value_in_unit(units.kilojoule_per_mole)

        delta_e_kjmol = final_e_kjmol - initial_e_kjmol
        msg = 'Minimized energy: {:8.3f} kJ/mol (deltaE = {:8.3} kJ/mol)'
        logging.info(msg.format(final_e_kjmol, delta_e_kjmol))

        if final_e_kjmol > 1000.0:  # arbitrary value
            emsg = 'Final energy is very high: {:8.3f} kJ/mol '
            emsg += '- minimization failed. Check structure for severe clashes.'
            raise StructureError(emsg.format(final_e_kjmol))

        self.__set_positions(state.getPositions())

        # Remove restraint force from system (LIFO)
        self._system.removeForce(self._system.getNumForces() - 1)
        self.potential_energy = final_e_kjmol

    # Neighbor Search
    # Thanks Xavier Martinez (IBPC, FR) for the KDTree suggestion.
    def get_neighbors(self, entity, radius=5.0, level='atom', method='exhaustive'):
        """Returns Atoms, Residues, or Chains in the vicinity of a given entity.

        Args:
            entity (:obj:): list or single instance of `Atom`, `Residue`, or
                `Chain` object(s).
            radius (float): distance threshold in Angstrom to consider an atom
                as neighbor of another. Default is 5.0 A.
            level (str): objects returned after search: 'atom', 'residue', or
                'chain'. Default is 'atom'.
            method (str): defines how to look for neighbors.
                Options are 'centroid' (calculates distance from center of
                gravity) or 'exhaustive' (looks for neighbors of all children of
                the object). Default is exhaustive.
        """

        if method not in ('centroid', 'exhaustive'):
            emsg = '\'method\' must be either: \'centroid\' or \'exhaustive\''
            raise ValueError(emsg)

        if level.lower() not in ('atom', 'residue', 'chain'):
            emsg = '\'level\' must be one of: \'atom\', \'residue\', or \'chain\''
            raise ValueError(emsg)
        else:
            level = level.lower()

        try:
            radius = float(radius)
        except ValueError as e:
            emsg = '\'radius\' should be a float, not {}'
            raise ValueError(emsg.format(type(radius))) from e

        if radius <= 0:
            raise ValueError('Distance threshold must be a positive number..')

        all_xyz = self._np_positions

        # Decompose into Atoms for search
        if isinstance(entity, app.topology.Atom):
            coords = [all_xyz[entity.index]]
            self_idx = {entity.index}
            method = 'exhaustive'

        elif isinstance(entity, (app.topology.Residue, app.topology.Chain)):
            if not hasattr(entity, 'atoms'):
                raise TypeError('Object must implement an \'atoms\' attribute.')
            atomidx = [a.index for a in entity.atoms()]
            coords = [all_xyz[idx] for idx in atomidx]
            self_idx = set(atomidx)

        elif isinstance(entity, list):  # list of any of the above?
            if not entity:  # empty list
                raise ValueError('You provided an empty list.')

            _types = (app.topology.Atom, app.topology.Residue, app.topology.Chain)
            list_types = sum([isinstance(item, _types) for item in entity])
            if list_types != len(entity):
                emsg = 'List of objects must contain Atoms, Residues, or Chains'
                raise TypeError(emsg)

            coords = []
            self_idx = set()
            for item in entity:
                if isinstance(item, app.topology.Atom):
                    coords.append(all_xyz[item.index])
                    self_idx.add(item.index)
                else:
                    atomidx = [a.index for a in item.atoms()]
                    coords += [all_xyz[idx] for idx in atomidx]
                    self_idx.update(atomidx)

            if len(coords) == 1:  # treat like single atom
                method = 'exhaustive'
                self_idx = {item.index}

        else:
            emsg = 'Object \'{}\' not supported in search. '
            emsg += 'Provide Atom(s), Residue(s), or Chain(s).'
            raise TypeError(emsg.format(entity))

        # Build KDTree if not there
        if self._kdt is None:
            self.__build_kdtree()

        # Perform search
        # kdt returns Point (p.index, p.radius) objects
        if method == 'centroid':
            _xyz = np.array(coords, dtype='d')  # redundant for Atom
            centroid = _xyz.mean(axis=0)
            assert centroid.shape == (3, ), \
                'Something went wrong when calculating object c.o.m.'

            neighbor_idx = set(a.index for a in self._kdt.search(centroid, radius))
            # neighbor_list = self._kdt.search(centroid, radius)
            # neighbor_idx = {item.index for item in neighbor_list}

        elif method == 'exhaustive':
            neighbor_idx = set()  # avoid duplicates
            for atom in coords:
                # results = (a.index for a in self._kdt.search(atom, radius))
                # neighbor_list = [it.index for it in results]
                neighbor_idx.update(a.index for a in self._kdt.search(atom, radius))

        # Exclude self neighbors
        # neighbor_idx = neighbor_idx.difference(self_idx)
        neighbor_idx -= self_idx  # difference update

        # Propagate to proper level (if not A)
        if level == 'atom':
            return [a for a in self.topology.atoms() if a.index in neighbor_idx]

        elif level == 'residue':
            result = []
            for residue in self.topology.residues():
                atomlist = {a.index for a in residue.atoms()}
                if atomlist & neighbor_idx:
                    result.append(residue)
            return result

        elif level == 'chain':
            result = []
            for chain in self.topology.chains():
                atomlist = {a.index for a in chain.atoms()}
                if atomlist & neighbor_idx:
                    result.append(chain)
            return result

    def get_neighboring_pairs(self, radius=5.0, level='atom'):
        """Returns all pairs of entities within a given radius of each other.
        """

        # Utility functions to retrieve parent objects
        def get_residue(atomdict, atom_idx):
            """Returns the Residue object associated with the atom index.

            Waives checking of index == None for performance. Assumes the
            output comes from the KDtree search so it should be fine.
            """
            atom = atomdict.get(atom_idx)
            return atom.residue

        def get_chain(atomdict, atom_idx):
            """Returns the Chain object associated with the atom index.

            Waives checking of index == None for performance. Assumes the
            output comes from the KDtree search so it should be fine.
            """
            atom = atomdict.get(atom_idx)
            return atom.residue.chain

        if level.lower() not in ('atom', 'residue', 'chain'):
            emsg = '\'level\' must \'atom\', \'residue\', or \'chain\''
            raise ValueError(emsg)
        else:
            level = level.lower()
            if level == 'residue':
                get_parent = get_residue
            elif level == 'chain':
                get_parent = get_chain

        try:
            radius = float(radius)
        except ValueError as e:
            vt = type(radius)
            emsg = '\'radius\' should be a float, not {}'
            raise ValueError(emsg.format(vt)) from e

        if radius <= 0:
            raise ValueError('Distance threshold must be a positive number..')

        msg = 'Searching all \'{}\' neighbor pairs within {} Angstrom'
        logging.debug(msg.format(level, radius))

        # Build KDTree if not there
        if self._kdt is None:
            self.__build_kdtree()

        raw_neighbors = self._kdt.neighbor_search(radius)

        # Filter according to requested level
        atomdict = {a.index: a for a in self.topology.atoms()}
        pairs_of_neighbors = []

        if level == 'atom':
            # Unpack structure to 3-item tuple
            unpacked = [(atomdict.get(n.index1),
                         atomdict.get(n.index2),
                         n.radius) for n in raw_neighbors]

        else:
            MAX_D = sys.maxsize
            min_distances = {}  # stores minimum distances between pairs
            for p in raw_neighbors:
                # Indices come sorted (p.index1 < p.index2)
                obj_i = get_parent(atomdict, p.index1)
                obj_j = get_parent(atomdict, p.index2)

                # exclude self-self
                if obj_i == obj_j:
                    continue

                obj_pair = (obj_i, obj_j)
                d_ij = p.radius

                cur_d = min_distances.get(obj_pair, MAX_D)
                if d_ij < cur_d:
                    min_distances[obj_pair] = d_ij

            unpacked = [(i, j, d) for (i, j), d in min_distances.items()]

        pairs_of_neighbors = unpacked

        # Return neighbors
        msg = 'Search returned {} pairs of neighbors'
        logging.debug(msg.format(len(pairs_of_neighbors)))
        return pairs_of_neighbors

    #
    # Energy-related functions
    #

    def add_energy(self, ene_term):
        """Method to couple an energy term to a `Structure` object

        Args:
            ene_term (:obj:`BaseEnergyTerm`): instance of `BaseEnergyTerm`.

        Raises:
            KeyError: if the energy term is not found.
            EnergyError: if energy term cannot be coupled to this `Structure`.
        """
        pass

    def remove_energy(self, eterm_name):
        """Method to remove a previously coupled energy term from a `Structure`.

        Args:
            eterm_name (str): name of the energy term, as defined in its class.

        Raises:
            KeyError: if the energy term is not coupled to the `Structure`.
        """
        pass
