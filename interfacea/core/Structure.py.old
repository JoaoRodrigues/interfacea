#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

Module containing Structure class to represent 3D molecular objects.
"""

import logging
import pathlib
import sys
import tempfile
import warnings

import numpy as np

import pdbfixer as pf
from pdbfixer.pdbfixer import Sequence

import simtk.openmm.app as app
import simtk.unit as units

from interfacea import RANDOM_SEED
from interfacea.chemistry import data
from interfacea.src.kdtree import kdtrees

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


# Classes

class StructureError(Exception):
    """Base class for all exceptions related to `Structure` objects."""
    pass


class StructureKDTreeError(StructureError):
    """Exception to signal errors in KDTree code/methods"""
    pass


class StructureWriteError(StructureError):
    """Exception to signal error when writing Structure to disk."""
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
        name (str): name of the Structure object.
        parameterized (bool): flags if the structure has been processed
            with a ForceField.
        positions (:obj:`OpenMM Positions`): OpenMM positions array.
        sequences (:obj:`list(`PDBFixer Sequence`)`: list of sequences
            contained in the `Structure` object.
        topology (:obj:`OpenMM Topology`): OpenMM topology.
    """

    def __init__(self, name, structure, build_kdtree=True):

        self.name = name
        self.parameterized = False

        self.topology = None
        self.positions = None

        # Cache objects
        self._kdt = None
        self._modeller = None
        self._pdbfixer = None

        # Internal/Private variables
        self._positions_np = None  # double-precision xyz np.array

        self.__set_positions(structure.positions)
        self.__set_topology(structure.topology)

        # Build KDTree
        if build_kdtree:
            self.__build_kdtree()

        logging.info('Created Structure from \'{}\''.format(name))

    def __repr__(self):
        """Print pretty things when called.
        """

        n = self.name
        nc = self.topology.getNumChains()
        nr = self.topology.getNumResidues()
        na = self.topology.getNumAtoms()

        return f"Structure '{n}' [{nc} chains | {nr} residues | {na} atoms]"

    def __add__(self, other):
        """Combines two Structure objects together"""

        this = app.Modeller(self.topology, self.positions)
        try:
            this.add(other.topology, other.positions)
        except Exception as e:
            emsg = f'Could not combined {this.name} and {other.name}'
            raise StructureError(emsg) from e

        return this

    # Pickling/Copying methods
    def copy(self):
        """Returns a (deep) copy of the Structure.

        Does not include KDTree or System objects to allow pickling.
        """

        if self.parameterized:
            warnings.warn('Structure copy is not parameterized')

        # KDTree is built on search anyway
        return self.__class__(self.name, self, build_kdtree=False)

    def __copy__(self, *args):
        """Shallow Copy Override.
        """

        return self.copy()

    def __deepcopy__(self, *args):
        """Deepcopy Override.
        """

        return self.copy()

    #
    # Private Methods
    #
    def __set_topology(self, topology):
        """Utility method to apply changes on topology changes.
        """

        logging.debug(f"Setting topology for '{self.name}'")
        self.topology = topology

    def __set_positions(self, positions):
        """Utility method to apply changes on atom addition/deletion.
        """

        logging.debug(f"Setting coordinates for '{self.name}'")
        self.positions = positions

        # Convert positions to numpy array
        # use double precision for coordinates because of KDTree
        _xyz_list = positions.value_in_unit(units.angstrom)
        self._positions_np = np.asarray(_xyz_list, dtype="d")
        del _xyz_list

        # Update KDTree
        if self._kdt is not None:
            self.__build_kdtree()

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
        try:
            self._kdt = kdtrees.KDTree(self._positions_np)
        except Exception as e:
            raise StructureKDTreeError from e

    def __load_to_pdbfixer(self):
        """Method to write a temporary PDB file and reload using PDBFixer.

        If PDBFixer was never called before, runs and caches the resulting
        Structure object. Always resets/empties the missing lists to avoid
        conflicts.
        """

        if self._pdbfixer is None:
            logging.debug('Creating PDBFixer object')

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

        self._pdbfixer.sequences = sequences
        self._pdbfixer.missingAtoms = {}
        self._pdbfixer.missingResidues = {}
        self._pdbfixer.missingTerminals = {}

        logging.debug('Cached PDBFixer Data Structure')

    def __load_openmm_modeller(self):
        """Method to create or load an OpenMM Modeller object.
        """

        if self._modeller is None:
            self._modeller = app.Modeller(self.topology, self.positions)

    #
    # Public Methods
    #

    def write(self, output, ftype=None, overwrite=False):
        """Writes structure to a file on disk.

        Uses OpenMM PDBFile or PDBxFile methods to write the `Structure` to a
        file on disk in PDB or mmCIF format, respectively. The output format is
        guessed from the user-provided file name or by the optional argument
        `format`.


        Args:
            output (str): name of file to write Structure to disk.
            ftype (str): file format to use when writing the structure to disk.
                Must be 'pdb' or 'cif'. If None (default), guesses the file
                type from the provided name.
            overwrite(bool, optional): force writing of the file, even if a
                file with the same name already exists. Default is False.

        Raises:
            StructureWriteError
        """

        _writers = {'cif': app.PDBxFile.writeFile,
                    'pdb': app.PDBFile.writeFile}

        try:
            fpath = pathlib.Path(output)
        except TypeError as e:
            emsg = f"'output' parameter should be a string: {type(output)}"
            raise StructureWriteError(emsg) from e

        # Guess file format
        if ftype is None:
            ftype = fpath.suffix.lstrip('.')
            if not ftype.strip():  # empty (no extension?)
                emsg = f'Failed to guess type from file name: {output}'
                raise StructureWriteError(emsg)
            logging.debug(f'Auto-assigned file type from file name: {ftype}')

        # Assign writer class
        writer = _writers.get(ftype)
        if writer is None:
            emsg = f"Unsupported or unknown file type: '{ftype}'"
            raise StructureWriteError(emsg)

        # Write
        if fpath.exists() and not overwrite:
            emsg = 'File already exists. Use overwrite=True or remove file.'
            raise StructureWriteError(emsg)

        try:
            with fpath.open(mode='w') as handle:
                writer(self.topology, self.positions, handle, keepIds=True)
        except Exception as e:
            emsg = f'Unexpected error when writing file: {output}'
            raise StructureWriteError(emsg) from e

        msg = f'Wrote structure to file: {fpath.resolve()}'
        logging.info(msg)

    # def delete(self, selection):
    #     """Removes atoms matching the selection from the Structure"""
    #     raise NotImplementedError

    def remove_waters(self, solvent_name=('HOH', 'WAT')):
        """Removes water molecules from the structure

        Arguments:
            solvent_name (tuple): three-letter residue names (as strings) to
                match internal Structure residue names. Default: HOH and WAT.
        """

        _names = set(solvent_name)

        self.__load_openmm_modeller()  # cache
        m = self._modeller

        matches = [r for r in m.topology.residues() if r.name in _names]
        m.delete(matches)

        self.__set_topology(m.topology)
        self.__set_positions(m.positions)

        msg = f'Removed {len(matches)} solvent molecules from Structure'
        logging.info(msg)

    def remove_unknowns(self, exceptions=None):
        """Removes all residues for which there is no template in the
        residue library.

        Arguments:
            exceptions (list): three-letter residue names for residues
                to skip even if they are not in the template library.
        """

        if exceptions is None:
            exceptions = []
        exceptions_set = set(exceptions)

        self.__load_to_pdbfixer()
        s = self._pdbfixer

        unknowns = []
        for r in self.topology.residues():
            if r.name not in s.templates and r.name not in exceptions_set:
                logging.debug(f'Deleting unknown residue: {r.name}')
                unknowns.append(r)

        if unknowns:
            self.__load_openmm_modeller()
            m = self._modeller
            m.delete(unknowns)

            self.__set_topology(m.topology)
            self.__set_positions(m.positions)

            msg = f'Removed {len(unknowns)} unknown residues from Structure'
            logging.info(msg)

    def add_residue_templates(self, template):
        """Adds templates to the internal residue library to allow
        building missing HEAVY atoms on non-standard residues.

        The template PDB File should have one residue only and should
        contain all bonding information in CONECT statements. The file
        should also be named after the residue you wish to add to the
        template library: e.g. ASP.pdb -> ASP residue.

        Hydrogen atoms are ignored. Parameterize your residue properly
        if you want to use them.

        Arguments:
            templates (str): path to PDB file to use as template to
                find/add missing heavy atoms.
        """

        # load cached pdbfixer object
        self.__load_to_pdbfixer()
        s = self._pdbfixer

        # Validate user input
        tpath = pathlib.Path(template)
        if not tpath.exists():
            emsg = f'Template PDB file not found: {template}'
            raise StructureError(emsg)

        tname = tpath.stem
        if tname in s.templates:
            wmsg = f'Residue {tname} is already in the template library'
            warnings.warn(wmsg)

        pdbfile = app.PDBFile(str(tpath.resolve()))

        # Remove Hs
        hydrogen = app.element.hydrogen
        _m = app.Modeller(pdbfile.topology, pdbfile.positions)
        hydrogens = [a for a in _m.topology.atoms() if a.element == hydrogen]
        _m.delete(hydrogens)

        pdbfile.topology = _m.topology
        pdbfile.positions = _m.positions

        s.templates[tname] = pdbfile
        logging.debug(f"Added '{tname}' to template library")

    def add_missing_heavy_atoms(self):
        """Adds heavy atoms that are missing from the structure.

        Relies on PDBFixer (topologies) to find and add missing atoms,
        namely the template library. If the Structure has residues not
        present in the PDBFixer residue library, add them from ideal
        coordinates using add_residue_templates().
        """

        self.__load_to_pdbfixer()
        s = self._pdbfixer

        s.findMissingResidues()
        s.findMissingAtoms()

        s.missingTerminals = {}  # do not add missing terminals here.

        num_missing = len(s.missingAtoms)

        if num_missing:
            logging.info(f'Adding {num_missing} missing heavy atoms')

            try:
                s.addMissingAtoms(seed=RANDOM_SEED)
            except TypeError:
                msg = ('The version of OpenMM you have installed does not add '
                       'missing atoms reproducibly - their positions may vary '
                       'There is nothing wrong with your structure. To remove '
                       'this warning, consider installing the latest version '
                       'of OpenMM from git')
                warnings.warn(msg)

                s.addMissingAtoms()

            self.__set_topology(s.topology)
            self.__set_positions(s.positions)

            # Issue warning about atom positions
            msg = ('Atoms added but their positions are not optimized. '
                   'Minimize the structure before doing any analysis')
            warnings.warn(msg)

    def add_capping_groups(self, chain_capping=None):
        """Adds capping groups to N- and C- termini of protein/nucleic acids.

        Use chain_capping to specify which capping groups you want to add to
        each individual chain in the Structure. By default, this method adds
        acetyl (ACE) and N-methyl (NME) groups to the N- and C- termini of
        proteins.

        Arguments:
            chain_capping (list): per-chain definition for capping groups. If
                not None, the number of items in the list must match the number
                of chains in the Structure. Allows options for each chain are
                'ACE', 'NME', or None (charged terminus).

        Raises:
            StructureError
        """

        _n_caps = set(('ACE', None))
        _c_caps = set(('NME', None))

        protein_aa = data.protein_aa

        self.__load_to_pdbfixer()
        s = self._pdbfixer

        chains = list(self.topology.chains())
        num_chains = len(chains)

        if chain_capping is None:
            chain_capping = [('ACE', 'NME') for _ in chains]
        else:
            # Validate number matches
            num_caps = len(chain_capping)
            if num_caps != num_chains:
                emsg = ("Number of capping groups must match number of chains:"
                        f" {num_caps} != {num_chains}")
                raise StructureError(emsg)

        # Find missing atoms that are NOT at the termini
        # to avoid adding them when calling this method.
        # We 'find' missing atoms now and store them, not to add them later.
        s.findMissingResidues()
        s.findMissingAtoms()
        missing_res = set(s.missingResidues.keys())
        missing_atm = set(s.missingAtoms.keys())

        # Build Sequence object(s) to add caps
        sequences = []
        for chain_idx, chain in enumerate(chains):

            # Validate capping choices
            n_cap, c_cap = chain_capping[chain_idx]
            if n_cap not in _n_caps:
                emsg = f"Unknown N-terminal capping group: {n_cap}"
                raise StructureError(emsg)

            if c_cap not in _c_caps:
                emsg = f"Unknown C-terminal capping group: {c_cap}"
                raise StructureError(emsg)

            # Get N- and C- terminal residues
            reslist = [r.name for r in chain.residues()]
            n_ter, c_ter = reslist[0], reslist[-1]

            # Do we need to add caps? Are the termini amino acids?
            if n_cap and n_ter != n_cap and n_ter in protein_aa:
                reslist.insert(0, n_cap)
                logging.debug(f"Capped chain {chain.id} with '{n_cap}' [N]")
            if c_cap and c_ter != c_cap and c_ter in protein_aa:
                reslist.append(c_cap)
                logging.debug(f"Capped chain {chain.id} with '{c_cap}' [C]")

            sequences.append(Sequence(chain.id, reslist))

        s.sequences = sequences  # set pdbfixer.sequences

        # Find missing heavy atoms in structure (after adding capping groups)
        # Remove previously found missing atoms
        s.findMissingResidues()
        for res in list(s.missingResidues.keys()):
            if res in missing_res:
                del s.missingResidues[res]

        s.findMissingAtoms()
        for res in list(s.missingAtoms.keys()):
            if res in missing_atm:
                del s.missingAtoms[res]

        # Add missing terminal atoms/residues
        try:
            s.addMissingAtoms(seed=RANDOM_SEED)
        except TypeError:
            msg = ('The version of OpenMM you have installed does not add '
                   'missing atoms reproducibly - their positions may vary '
                   'There is nothing wrong with your structure. To remove '
                   'this warning, consider installing the latest version '
                   'of OpenMM from git')
            warnings.warn(msg)

            s.addMissingAtoms()

        self.__set_topology(s.topology)
        self.__set_positions(s.positions)

    def mutate(self, residue, mutation):
        """Mutates a residue in the Structure to another.

        Works only for residues in the template library, and then likely
        only for protein/nucleic acid residues. The 'mutation' algorithm
        keeps atoms common to the original and mutated residue, so it
        preserves side-chain orientations if you are mutating amino acids.

        The positions of the added/replaced atoms are NOT optimized.

        Arguments:
            residue (:obj: Residue): Structure topology residue object
                to mutate.
            mutation (str): three-letter residue name for the mutated
                residue.

        Raises:
            StructureError
        """

        self.__load_to_pdbfixer()
        s = self._pdbfixer

        # Build list of valid residues to mutate to
        _pf = pf.pdbfixer
        resnames = _pf.proteinResidues + _pf.dnaResidues + _pf.rnaResidues
        _supported_resnames = set(resnames)

        if mutation not in _supported_resnames:
            emsg = f"Residue '{mutation}' not in template library"
            raise StructureError(emsg)

        # Do not add missing heavy atoms to other residues besides mutation.
        # applyMutations changes residue objects (uses Modeller) so we cannot
        # compare missingAtoms() directly.
        s.findMissingResidues()
        s.findMissingAtoms()
        missing_res = set(s.missingResidues.keys())  # tuple, no Res obj
        missing_atm = {(r.chain.id, r.id, r.name, r.insertionCode)
                       for r in s.missingAtoms.keys()}

        # Try to apply mutation
        # Defer to PDBFixer to catch other errors
        try:
            chain = residue.chain.id
            mut_list = [f'{residue.name}-{residue.id}-{mutation}']
            s.applyMutations(mut_list, chain)
        except Exception as e:
            emsg = f"Unknown error when applying mutation: {mut_list[0]}"
            raise StructureError(emsg) from e
        else:
            # Modify Structure
            s.findMissingResidues()
            for r in list(s.missingResidues.keys()):
                if r in missing_res:
                    del s.missingResidues[r]

            s.findMissingAtoms()
            for r in list(s.missingAtoms.keys()):
                r_uid = (r.chain.id, r.id, r.name, r.insertionCode)
                if r_uid in missing_atm:
                    del s.missingAtoms[r]

            s.missingTerminals = {}  # clear 'missing' termini

            try:
                s.addMissingAtoms(seed=RANDOM_SEED)
            except TypeError:
                msg = ('The version of OpenMM you have installed does not add '
                       'missing atoms reproducibly - their positions may vary '
                       'There is nothing wrong with your structure. To remove '
                       'this warning, consider installing the latest version '
                       'of OpenMM from git')
                warnings.warn(msg)

                s.addMissingAtoms()

        self.__set_topology(s.topology)
        self.__set_positions(s.positions)

        # Issue warning about atom positions
        wmsg = ('Residue mutated but atom positions were not optimized. '
                'Minimize the structure before doing any analysis')
        warnings.warn(wmsg)

    # Neighbor Search
    # Thanks Xavier Martinez (IBPC, FR) for the KDTree suggestion.
    def __unpack_entity(self, entity):
        """Returns a list of Atom indexes from (a list of Atom,) Residue,
        or Chain objects.

        Arguments:
            entity (:obj:): list or single instance of `Atom`, `Residue`,
                or `Chain` object(s).

        Returns:
            atom_idx_list (list): list containing the Atom indexes
                for each Atom in the entity.
        """

        _openmm_types = (app.topology.Residue, app.topology.Chain)

        if not isinstance(entity, list):
            _list = [entity]
        else:
            _list = entity

        # Unpack to Atoms
        atom_idx_list = []
        for e in _list:
            if isinstance(e, app.topology.Atom):
                atom_idx_list.append(e.index)
            elif isinstance(e, _openmm_types):
                for a in e.atoms():
                    atom_idx_list.append(a.index)
            else:
                emsg = f"Entity must be Residue, Chain: {type(e)}"
                raise ValueError(emsg)

        return atom_idx_list

    def __validate_neighbor_radius(self, radius):
        """Raises exceptions if the radius value is not a positive int"""

        try:
            radius = float(radius)
        except ValueError as e:
            emsg = f"radius must be a float, not {type(radius)}"
            raise ValueError(emsg) from e
        else:
            if radius <= 0.0:
                emsg = f"radius must be a positive number: {radius}"
                raise ValueError(emsg)

        return radius

    def get_neighbors(self, entity,
                      radius=5.0, level='atom', method='exhaustive'):
        """Returns Atoms, Residues, or Chains in the vicinity of a given entity.

        Arguments:
            entity (:obj:): list or single instance of `Atom`, `Residue`, or
                `Chain` object(s).
            radius (float): distance threshold in Angstrom to consider an atom
                as neighbor of another. Default is 5.0 A.
            level (str): objects returned after search: 'atom', 'residue', or
                'chain'. Default is 'atom'.
            method (str): defines how to look for neighbors.
                Options are 'centroid' (calculates distance from center of
                gravity) or 'exhaustive' (looks for neighbors of all children
                of the object). Default is exhaustive.
        """

        if method not in ('centroid', 'exhaustive'):
            emsg = f"'method' must be 'centroid' or 'exhaustive': {method}"
            raise ValueError(emsg)

        if level.lower() not in ('atom', 'residue', 'chain'):
            emsg = f"'level' must be 'atom', 'residue' or 'chain': {level}"
            raise ValueError(emsg)
        else:
            level = level.lower()
            logging.debug(f"neighbor search level={level}")

        radius = self.__validate_neighbor_radius(radius)

        # Decompose into Atoms for search
        if isinstance(entity, app.topology.Atom):
            atom_idx = [entity.index]

            if method != 'exhaustive':
                msg = f"'method' auto-set to 'exhaustive' for single Atom"
                logging.info(msg)
                method = 'exhaustive'

        else:
            atom_idx = self.__unpack_entity(entity)

        # Get xyz coordinates
        all_xyz = self._positions_np
        coords = [all_xyz[idx] for idx in atom_idx]
        self_idx = set(atom_idx)  # to avoid self contacts at the end

        # Build KDTree if not there
        if self._kdt is None:
            self.__build_kdtree()

        # Perform search
        # kdt returns Point (p.index, p.radius) objects
        if method == 'centroid':
            _xyz = np.array(coords, dtype='d')
            centroid = _xyz.mean(axis=0)

            try:
                result = self._kdt.search(centroid, radius)
            except Exception as e:
                emsg = f"Unexpected error in KDTree search for: {entity}"
                raise StructureKDTreeError(emsg) from e

            neighbor_idx = {a.index for a in result}

        elif method == 'exhaustive':
            neighbor_idx = set()  # avoid duplicates
            for atom in coords:
                try:
                    result = self._kdt.search(atom, radius)
                except Exception as e:
                    emsg = f"Unexpected error in KDTree search for: {entity}"
                    raise StructureKDTreeError(emsg) from e

                neighbor_idx.update({a.index for a in result})

        # Exclude self from neighbors
        neighbor_idx -= self_idx

        msg = f"Found {len(neighbor_idx)} atom neighbors within {radius}A"
        logging.debug(msg)

        # Propagate to proper level (if not A)
        atoms = {a for a in self.topology.atoms() if a.index in neighbor_idx}
        if level == 'atom':
            return sorted(atoms, key=lambda x: x.index)
        elif level == 'residue':
            residues = {a.residue for a in atoms}
            return sorted(residues, key=lambda x: (x.id, x.insertionCode))
        elif level == 'chain':
            chains = {a.residue.chain for a in atoms}
            return sorted(chains, key=lambda x: x.id)

    def get_neighbor_pairs(self, radius=5.0, level='atom'):
        """Returns all pairs of entities within a given radius of each other.

        Arguments:
            radius (float): distance threshold in Angstrom to consider an atom
                as neighbor of another. Default is 5.0 A.
            level (str): objects returned after search: 'atom', 'residue', or
                'chain'. Default is 'atom'.

        Returns:
            list of three-item tuples with all pairs of entities, type depends
            on level, within at least one atom within radius of each other.
            Tuples are of the form: (entity_i, entity_j, min distance).
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

        _get_parent = {'residue': get_residue,
                       'chain': get_chain}

        if level.lower() not in ('atom', 'residue', 'chain'):
            emsg = f"'level' must be 'atom', 'residue' or 'chain': {level}"
            raise ValueError(emsg)
        else:
            level = level.lower()
            logging.debug(f"neighbor search level={level}")

        radius = self.__validate_neighbor_radius(radius)

        # Build KDTree if not there
        if self._kdt is None:
            self.__build_kdtree()

        try:
            raw_neighbors = self._kdt.neighbor_search(radius)
        except Exception as e:
            emsg = f"Unexpected error in KDTree neighbor search"
            raise StructureKDTreeError(emsg) from e

        # Filter according to requested level
        atomdict = {a.index: a for a in self.topology.atoms()}
        result = []

        if level == 'atom':
            # Unpack structure to 3-item tuple
            result = [(atomdict.get(n.index1),
                       atomdict.get(n.index2),
                       n.radius) for n in raw_neighbors]

        else:
            get_parent = _get_parent.get(level)

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

            result = [(i, j, d) for (i, j), d in min_distances.items()]

        # Return neighbors
        msg = f"Found {len(result)} neighbor pairs within {radius}A"
        logging.debug(msg)

        return result
