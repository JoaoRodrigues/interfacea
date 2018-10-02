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

import networkx as nx

import pdbfixer as pf
from pdbfixer.pdbfixer import Sequence

import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as units

# Setup logger
# _private name to prevent collision/confusion with parent logger
_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())

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

    Attributes:
        topology (:obj:`OpenMM Topology`): OpenMM topology.
        positions (:obj:`OpenMM Positions`): OpenMM positions array.

        sequences(:obj:`list(`PDBFixer Sequence`)`: list of sequences described in
            the `Structure` object.
        forcefield(:obj:`ForceField`): pointer to associated `ForceField` class
            defined at runtime.
    """

    def __init__(self, name, structure):

        self.name = name
        self.topology = structure.topology
        self.positions = structure.positions

        self.sequences = None
        self.forcefield = None  # forcefield name (str)
        self._forcefield = None  # forcefield object

        self._system = None

        self._pdbfixer = None  # cache PDBFixer structure if we need it.

        self._bonds_per_atom = False  # flag set when _bonded_atoms is called
        self._residue_graphs = False  # flag set when _make_res_graphs is called

        _logger.debug('Created Structure from \'{}\''.format(name))

    def __repr__(self):
        """Print pretty things when called.
        """

        rep_str = 'Structure created from {}'.format(os.path.basename(self.name))

        n_chain = self.topology.getNumChains()
        n_resid = self.topology.getNumResidues()
        n_atoms = self.topology.getNumAtoms()

        rep_str += ' with {} chain, {} residues, and {} atoms'.format(n_chain, n_resid, n_atoms)

        if self._forcefield:
            rep_str += ' (ff={})'.format(self.forcefield)

        return rep_str

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
        _logger.debug('Adding termini to {} chains'.format(num_chains))

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
        # export them properly. For now, we save as PDB and re-read.
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
                _logger.debug('Capping chain {} N-ter with \'{}\''.format(chain.id, n_cap))
            if c_cap and c_ter != c_cap:
                chain_reslist.append(c_cap)
                _logger.debug('Capping chain {} C-ter with \'{}\''.format(chain.id, c_cap))

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
        s.addMissingAtoms()

        self.topology = s.topology
        self.positions = s.positions

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
        _logger.debug('Adding {} missing atoms'.format(n_added_atoms))

        if n_added_atoms:
            s.addMissingAtoms()
            self.topology = s.topology
            self.positions = s.positions

            # Issue warning about atom positions
            warnings.warn(('Atoms added but their positions are not optimized. '
                           'Make sure to minimize the structure before doing any analysis'))

    def protonate(self, forcefield='amber14-all.xml', pH=7.0):
        """Method to add hydrogen atoms to a `Structure` object.

        Uses the `Modeller` class from OpenMM to add hydrogen atoms to the structure.
        Removes existing hydrogen atoms (to avoid naming issues) before adding them
        again with naming and topology matching the force field and chosen pH.


        Args:
            forcefield (str): name of file defining the force field.
            pH (float): numerical value of the pH to check protonation states
                of ionizable groups.
        """

        if self._forcefield is None:
            self._load_forcefield(forcefield)

        model = app.Modeller(self.topology, self.positions)

        # Remove existing hydrogens first, just in case
        _elem_H = app.element.hydrogen
        existing_H = [a for a in model.topology.atoms() if a.element == _elem_H]
        model.delete(existing_H)

        _logger.debug('Protonating structure at pH {}'.format(pH))
        model.addHydrogens(forcefield=self._forcefield, pH=pH)

        self.topology = model.topology
        self.positions = model.positions

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

            mut_per_chain[chain].append('{}-{}-{}'.format(name, idx, new))

        # Mutate on each chain at a time
        for chain in mut_per_chain:
            muts = mut_per_chain[chain]
            _logger.debug('Making {} mutations on chain {}'.format(len(muts), chain))

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

            s.addMissingAtoms()

        self.topology = s.topology
        self.positions = s.positions

        # Issue warning about atom positions
        warnings.warn(('Residue mutated but atom positions are not optimized. '
                       'Make sure to minimize the structure before doing any analysis'))

    def write(self, output, ftype='cif', overwrite=False):
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
            raise TypeError('\'output\' argument must be \'str\' or a file-like object')

        try:
            with handle:
                writer(self.topology, self.positions, handle, keepIds=True)
        except Exception as e:
            emsg = 'Error when writing Structure to file: {}'.format(handle.name)
            raise StructureError(emsg) from e

    #
    # MM Functions
    #

    def _load_forcefield(self, forcefield='amber14-all.xml'):
        """Utility private method to load forcefield definitions.
        """

        try:
            loaded_forcefield = app.ForceField(forcefield)
            _logger.debug('Loaded forcefield definitions from: {}'.format(forcefield))

        except ValueError as e:
            emsg = 'Error when loading forcefield XML file: {}'.format(forcefield)
            raise StructureError(emsg) from e

        self._forcefield = loaded_forcefield
        self.forcefield = forcefield

    def parameterize(self, forcefield='amber14-all.xml'):
        """Wrapper function to create an OpenMM system from the Structure.

        By default, does not use a cutoff (nor I plan to change this).

        Args:
            forcefield (str): name of 'xml' file containing forcefield definitions.
        """

        if self._forcefield is None:
            self._load_forcefield(forcefield)

        system = self._forcefield.createSystem(self.topology, nonbondedMethod=app.NoCutoff)
        self._system = system
        _logger.debug('Successfully created new OpenMM System for structure')

    def _get_bonded_atoms(self):
        """Creates a dictionary of bonds per atom.
        """

        bond_dict = {}
        for residue in self.topology.residues():
            for b in residue.bonds():

                bond_dict.setdefault(b.atom1, [])
                bond_dict.setdefault(b.atom2, [])

                bond_dict[b.atom1].append(b.atom2)
                bond_dict[b.atom2].append(b.atom1)

            residue.bonds_per_atom = bond_dict

        _logger.debug('Created per-atom bonding dictionary from topology')
        self._bonds_per_atom = True

    def _make_residue_graphs(self):
        """Creates a networkx.Graph representation of a `Residue`.

        Uses the atom elements as node attributes and the topology bonds
        as edges.
        """

        for residue in self.topology.residues():
            at_to_idx = {at: idx for idx, at in enumerate(residue.atoms())}

            # Make graph of residue
            res_g = nx.Graph()
            for atom, idx in at_to_idx.items():
                res_g.add_node(idx, element=atom.element.atomic_number)

            for b in residue.bonds():
                a1, a2 = b.atom1, b.atom2
                a1_idx, a2_idx = at_to_idx.get(a1), at_to_idx.get(a2)
                if a1_idx is not None and a2_idx is not None:
                    res_g.add_edge(a1_idx, a2_idx)

            residue._g = res_g

        _logger.debug('Converted residue topologies to graph representation')
        self._residue_graphs = True

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