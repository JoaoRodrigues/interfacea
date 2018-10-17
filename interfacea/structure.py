#!/usr/bin/env python

"""
Analysis of Biomolecular Interfaces.

Module containing parsing and setup routines for
molecular structures.
"""

from __future__ import print_function

import logging
import io
import os
import tempfile
import warnings

import numpy as np
import networkx as nx

import pdbfixer as pf
from pdbfixer.pdbfixer import Sequence

import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as units

from src import kdtrees

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())

# .minimize([posre=True, nsteps=50])


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
        random_seed (int): integer to use as a random number seed.
            Default is 917.
        build_kdtree (bool): automatically build KDTree on instantiation.
            Default is True.

    Attributes:
        topology (:obj:`OpenMM Topology`): OpenMM topology.
        positions (:obj:`OpenMM Positions`): OpenMM positions array.

        sequences(:obj:`list(`PDBFixer Sequence`)`: list of sequences described in
            the `Structure` object.
        forcefield(:obj:`ForceField`): pointer to associated `ForceField` class
            defined at runtime.

        potential_energy (float): potential energy (in kJ/mol) calculated using
            the forcefield parameters.
    """

    def __init__(self, name, structure, random_seed=917, build_kdtree=True):

        self.sequences = None
        self.forcefield = None  # forcefield name (str)
        self.potential_energy = None

        self._forcefield = None  # forcefield object
        self._system = None
        self._pdbfixer = None  # cache PDBFixer structure if we need it.
        self._kdt = None

        self.seed = random_seed  # to allow reproducibility of results across the library

        self.name = name
        self._set_topology(structure.topology)
        self._set_positions(structure.positions)

        # Build KDTree
        if build_kdtree:
            self._build_kdtree()

        logging.debug('Created Structure from \'{}\''.format(name))

    def __repr__(self):
        """Print pretty things when called.
        """

        rep_str = 'Structure ({})'.format(os.path.basename(self.name))

        n_chain = self.topology.getNumChains()
        n_resid = self.topology.getNumResidues()
        n_atoms = self.topology.getNumAtoms()

        rep_str += ' ({} chain(s), {} residue(s), and {} atom(s))'.format(n_chain, n_resid, n_atoms)

        if self._forcefield:
            rep_str += ' (ff={})'.format(self.forcefield)

        return rep_str

    # Pickling/Copying methods
    def copy(self):
        """Returns a (deep) copy of the Structure.
        """

        newstruct = self.__class__(self.name, self, build_kdtree=False)
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

    def _build_kdtree(self):
        """Builds a KDTree for fast neighbor lookup.

        KDTree in C and Python bindings by Michiel de Hoon, taken from Biopython.
        For details, read the source code at src/kdtrees.c

        The module implements a KDTree class that takes a Nx3 numpy array of doubles.
        The class exposes three methods:
            search: returns indices of atoms within a radius (in Angstrom) around a central point.
            neighbor_search: returns all point pairs within a certain radius of each other.
            neighbor_search_simple: same as above, slow implementation for test purposes.

        In interfacea, we will expose only search but wrap it to return either Atom, Residue, or Chain objects,
        a bit like Biopython does it.
        """

        logging.debug('Building KDTree (this might take a minute or two)')
        self._kdt = kdtrees.KDTree(self._np_positions)

    def _set_topology(self, topology):
        """Utility method to apply changes on topology changes.
        """

        self.topology = topology

        self._get_bonded_atoms()
        self._make_residue_graphs()

    def _set_positions(self, positions):
        """Utility method to apply changes on atom addition/deletion.
        """

        self.positions = positions

        # Convert positions to numpy array
        _xyz_list = positions.value_in_unit(units.angstrom)
        self._np_positions = np.asarray(_xyz_list, dtype="d")  # double precision needed for kdtree
        del _xyz_list

        if self._kdt is not None:
            self._build_kdtree()

    def _load_to_pdbfixer(self):
        """Utility class to write a temporary PDB file and reload using PDBFixer.

        If PDBFixer was never called before, runs and caches the resulting Structure.
        Always resets/empties the missing lists to avoid conflicts.
        """

        if self._pdbfixer is None:
            with tempfile.TemporaryFile(mode='r+') as handle:
                app.PDBFile.writeFile(self.topology, self.positions, handle, keepIds=True)
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

        logging.debug('Cached PDBFixer data structure')

    def _get_bonded_atoms(self):
        """Creates a dictionary of bonds per atom.

        Somewhat a performance bottleneck for large structures.
        Should we optimize? How?
        """

        reslist = list(self.topology.residues())
        for residue in reslist:
            bond_dict = {}
            for b in residue.bonds():

                bond_dict.setdefault(b.atom1, [])
                bond_dict.setdefault(b.atom2, [])

                bond_dict[b.atom1].append(b.atom2)
                bond_dict[b.atom2].append(b.atom1)

            residue.bonds_per_atom = bond_dict

        logging.debug('Created per-atom bonding dictionary from topology')

    def _make_residue_graphs(self):
        """Creates a networkx.Graph representation of a `Residue`.

        Uses the atom elements as node attributes and the topology bonds
        as edges.

        Somewhat a performance bottleneck for large structures.
        Should we optimize? How?
        """

        reslist = list(self.topology.residues())
        for residue in reslist:
            at_to_idx = {at: idx for idx, at in enumerate(residue.atoms())}

            # Make graph of residue
            res_g = nx.Graph()
            for atom, idx in at_to_idx.items():
                res_g.add_node(idx, element=atom.element.atomic_number)
                for bonded in residue.bonds_per_atom[atom]:
                    if bonded in at_to_idx:
                        res_g.add_edge(idx, at_to_idx[bonded])

            residue._g = res_g

        logging.debug('Converted residue topologies to graph representation')

    def _load_forcefield(self, forcefield='amber14-all.xml'):
        """Utility private method to load forcefield definitions.
        """

        try:
            loaded_forcefield = app.ForceField(forcefield)
            logging.debug('Loaded forcefield definitions from: {}'.format(forcefield))

        except ValueError as e:
            emsg = 'Error when loading forcefield XML file: {}'.format(forcefield)
            raise StructureError(emsg) from e

        self._forcefield = loaded_forcefield
        self.forcefield = forcefield

    #
    # Public Methods
    #

    # IO
    def write(self, output, ftype=None, overwrite=False):
        """Writes `Structure` object to file.

        Uses OpenMM PDBFile or PDBxFile methods to write the `Structure` to a file on disk in
        PDB or mmCIF format, respectively. The output format is guessed from the user-provided
        file name or by the optional argument `format`.


        Args:
            output (file/str): file object or name to create the new file on disk.
            ftype (str): format to use when writing the file. Must be either 'pdb' or 'cif'.
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
                emsg = 'You must either provide an extension or a filename with one.'
                raise StructureError(emsg)

        writer = _writers.get(ftype)
        if writer is None:
            emsg = 'Unsupported file type \'{}\'. Choose from {}'.format(ftype, _fmt_str)
            raise StructureError(emsg)

        if isinstance(output, str):
            if os.path.isfile(output) and not overwrite:
                emsg = 'File already exists. Use overwrite=True or remove file.'
                raise OSError(emsg)
            handle = open(output, 'w')

        elif isinstance(output, io.IOBase) or (hasattr(output, 'file') and
                                               isinstance(output.file, io.IOBase)):
            handle = output
        else:
            raise TypeError('\'output\' argument must be a file name or a file-like object')

        try:
            with handle:
                writer(self.topology, self.positions, handle, keepIds=True)
        except Exception as e:
            emsg = 'Error when writing Structure to file: {}'.format(handle.name)
            raise StructureError(emsg) from e

    # Structure manipulation
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
        """Method to add missing terminal atoms to (all) chains in the `Structure`.

        Uses PDBFixer to modify structure. By default, adds acetyl (ACE) and N-methyl (NME) caps
        to N- and C- termini, respectively. Other caps can be specified using the `ends` option.

        Args:
            ends (:obj:`list(tuple(str, str)`, optional): definitions of termini groups to add
                to each chain. Number of items in list must match the number of chains in the
                `Structure`. Allowed options for each chain are: 'ACE', 'NME', or None (charged
                terminus).

        Raises:
            StructureError
        """

        _allowed_n_caps = set(('ACE', None))
        _allowed_c_caps = set(('NME', None))

        chains = list(self.topology.chains())
        num_chains = len(chains)
        logging.debug('Adding termini to {} chains'.format(num_chains))

        if ends is not None:
            num_ends = len(ends)
            if num_ends != num_chains:
                emsg = 'Number of terminal capping groups ({}) != '.format(num_ends)
                emsg += 'number of chains in the molecule ({})'.format(num_chains)
                raise StructureError(emsg)

            for chain_idx, chain in enumerate(ends):
                name = chains[chain_idx].id
                n_cap, c_cap = chain

                if n_cap not in _allowed_n_caps:
                    emsg = 'User-specified N-terminal for chain {}: {}'.format(name, n_cap)
                    raise StructureError(emsg)

                if c_cap not in _allowed_c_caps:
                    emsg = 'User-specified C-terminal for chain {}: {}'.format(name, c_cap)
                    raise StructureError(emsg)

        else:
            ends = [('ACE', 'NME') for _ in range(num_chains)]

        # PDBFixer is picky with mmCIF files and OpenMM does not
        # export them properly. The irony I know.
        # For now, we save as PDB and re-read.
        self._load_to_pdbfixer()
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

            # Add caps if necessary
            n_cap, c_cap = ends[chain_idx]
            if n_cap and n_ter != n_cap:
                chain_reslist.insert(0, n_cap)
                logging.info('Adding \'{}\' capping group to chain {} N-terminus'.format(n_cap, chain.id))
            if c_cap and c_ter != c_cap:
                chain_reslist.append(c_cap)
                logging.info('Adding \'{}\' capping group to chain {} C-terminus'.format(c_cap, chain.id))

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
        s.addMissingAtoms(seed=self.seed)

        self._set_topology(s.topology)
        self._set_positions(s.positions)

    def add_missing_atoms(self):
        """Method to add missing atoms to a `Structure` object.

        Uses PDBFixer to analyze structure and add missing heavy atoms.
        """

        # Save as PDB and re-read.
        # Someone should really write an OpenMM to PDBFixer conversion ...
        self._load_to_pdbfixer()
        s = self._pdbfixer

        s.findMissingResidues()
        s.findMissingAtoms()

        s.missingTerminals = {}  # do not add missing terminals here.

        n_added_atoms = len(s.missingAtoms)

        if n_added_atoms:
            logging.info('Found missing heavy atoms: {}'.format(n_added_atoms))
            s.addMissingAtoms(seed=self.seed)
            self._set_topology(s.topology)
            self._set_positions(s.positions)

            # Issue warning about atom positions
            warnings.warn(('Atoms added but their positions are not optimized. '
                           'Make sure to minimize the structure before doing any analysis'))

    def protonate(self, forcefield='amber14-all.xml', pH=7.0, keep_existing=False):
        """Method to add hydrogen atoms to a `Structure` object.

        Uses the `Modeller` class from OpenMM to add hydrogen atoms to the structure.
        Removes existing hydrogen atoms (to avoid naming issues) before adding them
        again with naming and topology matching the force field and chosen pH.


        Args:
            forcefield (str): name of file defining the force field.
            pH (float): numerical value of the pH to check protonation states
                of ionizable groups.
            keep_existing (bool): Default is False. Does not remove existing protons.
                This is dangerous because protons might not follow proper naming conventions
                but useful in the case of wanting specific protonation states.
        """

        if self._forcefield is None:
            self._load_forcefield(forcefield)

        model = app.Modeller(self.topology, self.positions)

        if not keep_existing:
            _elem_H = app.element.hydrogen
            existing_H = [a for a in model.topology.atoms() if a.element == _elem_H]
            model.delete(existing_H)

        logging.info('Protonating structure at pH {}'.format(pH))
        model.addHydrogens(forcefield=self._forcefield, pH=pH)

        self._set_topology(model.topology)
        self._set_positions(model.positions)

        # Issue warning about atom positions
        warnings.warn(('Protons added but their positions are not optimized. '
                       'Make sure to minimize the structure before doing any analysis'))

    def mutate(self, mutation_list):
        """Mutates residues in the molecule using PDBFixer.

        This is a very crude method of deleting/adding atoms, so most useful (or reasonable)
        for single mutations. Mutations of multiple residues at once might yield a
        very bad structure, even if followed by minimization.

        Args:
            mutation_list (:obj:`list(tuple)`): list of two-item tuples containing the
                id of the residue to mutate as a string with 'chain-resname-resid' and
                the resname of the mutated residue. Residue names should always be in
                three-letter code to avoid ambiguities.
                e.g. ('A-ASN-1', 'ALA') mutates ASN1 of chain A to alanine.

        Raises:
            StructureError
        """

        self._load_to_pdbfixer()
        s = self._pdbfixer

        # Mutate only. Do not add/complete structure.
        s.findMissingResidues()
        s.findMissingAtoms()
        missing_residues = set(s.missingResidues.keys())
        missing_atoms = set(s.missingAtoms.keys())

        # Build list of valid residues to mutate to
        _supported_resnames = set(pf.pdbfixer.proteinResidues +
                                  pf.pdbfixer.dnaResidues +
                                  pf.pdbfixer.rnaResidues)

        # Sanity check on mutation list
        mut_per_chain = {}
        for mutation in mutation_list:
            try:
                ori, new = mutation
                chain, name, idx = ori.split('-')
            except Exception as e:
                emsg = 'Wrong format in mutation: \'{}\''.format(mutation)
                raise StructureError(emsg) from e

            if new not in _supported_resnames:
                emsg = 'Residue not supported for mutation: {}'.format(new)
                raise StructureError(emsg)

            # defer to PDBFixer to catch errors of mutating non-existing residues
            # or on non-existing chains.

            if chain not in mut_per_chain:
                mut_per_chain[chain] = []

            logging.info('Mutating [{}]{}{} to {}'.format(name, idx, chain, new))
            mut_per_chain[chain].append('{}-{}-{}'.format(name, idx, new))

        # Mutate on each chain at a time
        for chain in mut_per_chain:
            muts = mut_per_chain[chain]

            try:
                s.applyMutations(muts, chain)
            except (KeyError, ValueError) as e:
                emsg = 'There was an error when applying mutations to the structure'
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

            s.addMissingAtoms(seed=self.seed)

        self._set_topology(s.topology)
        self._set_positions(s.positions)

        # Issue warning about atom positions
        warnings.warn(('Residue mutated but atom positions are not optimized. '
                       'Make sure to minimize the structure before doing any analysis'))

    # MM Functions

    def parameterize(self, forcefield='amber14-all.xml'):
        """Wrapper function to create an OpenMM system from the Structure.

        By default, does not use a cutoff (nor I plan to change this).

        Args:
            forcefield (str): name of 'xml' file containing forcefield definitions.
        """

        if self._forcefield is None or forcefield != self.forcefield:
            if forcefield != self.forcefield:
                warnings.warn(('Structure previously parameterized with \'{}\'. '
                               'It is advisable to run protonate() again'.format(self.forcefield)))

            self._load_forcefield(forcefield)

        system = self._forcefield.createSystem(self.topology, nonbondedMethod=app.NoCutoff)
        self._system = system

        logging.info('Structure parameterized using \'{}\''.format(forcefield))

    def calculate_energy(self):
        """Calculates the potential energy of the system.
        """

        if self._system is None:
            emsg = 'Structure is not parameterized yet. Use `parameterize()`.'
            raise StructureError(emsg)

        # Set integrator
        integrator = mm.LangevinIntegrator(300 * units.kelvin, 1.0 / units.picosecond,
                                           2.0 * units.femtosecond)
        integrator.setRandomNumberSeed(self.seed)
        integrator.setConstraintTolerance(0.00001)

        # Create context
        context = mm.Context(self._system, integrator)
        context.setPositions(self.positions)

        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(units.kilojoule_per_mole)
        logging.info('Potential energy of the structure = {:8.3f} kJ/mol'.format(energy))

        if energy > 0.0:  # sort of arbitrary value
            warnings.warn(('Potential energy of the structure is high ({:8.3f}). '
                           'You should try minimize() before any analysis'.format(energy)))

        self.potential_energy = energy

        del context
        del integrator

    def minimize(self, iterations=250, hydrogen_only=False):
        """Perform energy minimization using OpenMM.

        Args:
            iterations (int): number of steps taken by the minimizer (or 0 until convergence)
            hydrogen_only (bool): minimize positions of hydrogen atoms only. Default is False.
        """

        if self._system is None:
            emsg = 'Structure is not parameterized yet. Use `parameterize()`.'
            raise StructureError(emsg)

        # Set integrator
        integrator = mm.LangevinIntegrator(300 * units.kelvin, 1.0 / units.picosecond,
                                           2.0 * units.femtosecond)
        integrator.setRandomNumberSeed(self.seed)
        integrator.setConstraintTolerance(0.00001)

        # Harmonic position restraints
        posre = mm.CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
        posre.addGlobalParameter("k", 500.0 * (units.kilojoule_per_mole / units.nanometer**2))
        posre.addPerParticleParameter("x0")
        posre.addPerParticleParameter("y0")
        posre.addPerParticleParameter("z0")

        if hydrogen_only:

            logging.info('Optimizing positions of hydrogen atoms only')

            elemlist = [a.element.atomic_number for a in self.topology.atoms()]

            for idx, atom_crd in enumerate(self.positions):
                elem = elemlist[idx]
                if elem != 1:
                    posre.addParticle(idx, atom_crd.value_in_unit(units.nanometers))

            self._system.addForce(posre)

            n_restraints = posre.getNumParticles()
            n_atoms = len(elemlist)
            logging.debug('Restrained {} out of {} atoms'.format(n_restraints, n_atoms))

        # Perform minimization
        context = mm.Context(self._system, integrator)
        context.setPositions(self.positions)

        state = context.getState(getEnergy=True)
        initial_e = state.getPotentialEnergy().value_in_unit(units.kilojoule_per_mole)
        logging.debug('Energy before minimization: {:8.3f} kJ/mol'.format(initial_e))

        if initial_e > 1000000.0:  # sort of arbitrary value
            warnings.warn(('Potential energy of the starting structure is very high ({:8.3f} kJ/mol). '
                           'Minimization is likely to fail.'.format(initial_e)))

        mm.LocalEnergyMinimizer.minimize(context, maxIterations=iterations)

        state = context.getState(getPositions=True, getEnergy=True)
        final_e = state.getPotentialEnergy().value_in_unit(units.kilojoule_per_mole)

        delta_e = final_e - initial_e
        logging.info('Energy after minimization: {:8.3f} kJ/mol (deltaE = {:8.3} kJ/mol)'.format(final_e, delta_e))

        if final_e > 1000.0:  # arbitrary value
            emsg = 'Minimization failed. Minimized energy is still very high: {:8.3f} kJ/mol'
            raise StructureError(emsg.format(final_e))

        self._set_positions(state.getPositions())

        # Remove restraint force from system (LIFO)
        self._system.removeForce(self._system.getNumForces() - 1)

        self.potential_energy = final_e

    # Neighbor Search
    def get_neighbors(self, entity, radius=5.0, level='atom', method='exhaustive'):
        """Returns [A]toms, [R]esidues, or [C]hains in the vicinity of a given entity.

        Args:
            entity (:obj:): list or single instance of `Atom`, `Residue`, or `Chain` object(s).
            radius (float): distance threshold in Angstrom to consider an atom as neighbor.
                Default is 5.0 A.
            level (str): type of object returned after search: 'atom', 'residue', or 'chain'.
                Default is 'Atom'.
            method (str): defines how to look for neighbors.
                Options are 'centroid' (calculates distance from center of gravity) or
                'exhaustive' (looks for neighbors of all children of the object).
                Default is exhaustive.
        """

        msg = 'Search Parameters: radius={}/level={}/method={}'.format(radius, level, method)
        logging.debug(msg)

        if method not in ('centroid', 'exhaustive'):
            emsg = '\'method\' argument must be either: \'centroid\' or \'exhaustive\''
            raise ValueError(emsg)

        if level.lower() not in ('atom', 'residue', 'chain'):
            emsg = '\'level\' argument must be one of: \'atom\', \'residue\', or \'chain\''
        else:
            level = level.lower()

        try:
            radius = float(radius)
        except ValueError as e:
            raise ValueError('\'radius\' should be a float, not {}'.format(type(radius))) from e

        if radius <= 0:
            raise ValueError('Distance threshold must be a positive number ...')

        # Decompose into Atoms for search
        if isinstance(entity, app.topology.Atom):
            coords = [self._np_positions[entity.index]]
            self_idx = {entity.index}
            method = 'exhaustive'

        elif isinstance(entity, (app.topology.Residue, app.topology.Chain)):
            if not hasattr(entity, 'atoms'):
                raise TypeError('Object should implement an \'atoms\' attribute.')
            coords = [self._np_positions[a.index] for a in entity.atoms()]
            self_idx = {a.index for a in entity.atoms()}

        elif isinstance(entity, list):  # list of any of the above?
            if not entity:  # empty list
                raise ValueError('You provided an empty list.')

            _types = (app.topology.Atom, app.topology.Residue, app.topology.Chain)
            list_types = sum([isinstance(item, _types) for item in entity])  # True == 1
            if list_types != len(entity):
                raise TypeError('List of objects must contain only Atoms, Residues, or Chains')

            coords = []
            self_idx = set()
            for item in entity:
                if isinstance(item, app.topology.Atom):
                    coords.append(self._np_positions[item.index])
                    self_idx.add(item.index)
                else:
                    coords += [self._np_positions[a.index] for a in item.atoms()]
                    self_idx.update((a.index for a in item.atoms()))

            if len(coords) == 1:  # treat like single atom
                method = 'exhaustive'
                self_idx = {item.index}

        else:
            emsg = 'Object \'{}\' not supported in search. Provide Atom(s), Residue(s), or Chain(s).'.format(entity)
            raise TypeError(emsg)

        logging.debug('Search object comprises {} atoms'.format(len(coords)))
        # Perform search
        if method == 'centroid':
            logging.debug('Using \'centroid\' method to search')
            _xyz = np.array(coords, dtype='d')  # redundant for Atom but we pass...
            centroid = _xyz.mean(axis=0)
            assert centroid.shape == (3, ), 'Something went wrong when calculating object c.o.m.'

            neighbor_list = self._kdt.search(centroid, radius)  # list of Point (p.index, p.radius) objects
            neighbor_idx = {item.index for item in neighbor_list}

        elif method == 'exhaustive':
            logging.debug('Using \'exhaustive\' method to search')
            neighbor_idx = set()  # avoid duplicates
            for atom in coords:
                neighbor_list = [it.index for it in self._kdt.search(atom, radius)]
                neighbor_idx.update(neighbor_list)

        # Exclude self neighbors
        neighbor_idx = neighbor_idx.difference(self_idx)
        logging.debug('Search returned {} neighboring atoms'.format(len(neighbor_idx)))

        # Propagate to proper level (if not A)
        if level == 'atom':
            logging.debug('Returning Atom objects')
            return [a for a in self.topology.atoms() if a.index in neighbor_idx]

        elif level == 'residue':
            logging.debug('Returning Residue objects')
            result = []
            for residue in self.topology.residues():
                atomlist = {a.index for a in residue.atoms()}
                if atomlist & neighbor_idx:
                    result.append(residue)
            return result

        elif level == 'chain':
            logging.debug('Returning Chain objects')
            result = []
            for chain in self.topology.chains():
                atomlist = {a.index for a in chain.atoms()}
                if atomlist & neighbor_idx:
                    result.append(chain)
            return result

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