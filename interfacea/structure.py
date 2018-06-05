#!/usr/bin/env python

"""
Energy Evaluation of Biomolecular Interfaces.

Module containing parsing and setup routines for
molecular structures.
"""

from __future__ import print_function

import copy
import logging
import os
import tempfile

import pdbfixer as pf
from pdbfixer.pdbfixer import Sequence

import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as units

# Write mmCIF
# Write PDB

# keepIDs

# .add_termini()
# .add_missing_atoms()
# .mutate()

# .forcefield()
# .minimize([posre=True, nsteps=50])


class StructureError(Exception):
    """Dummy catch-all class for all exceptions related to `Structure` objects.

    Allows re-raising of exceptions, since exception chaining is not supported
    in Python 2.x.
    """

    def __init__(self, message, cause=None):
        if cause is not None:
            super(StructureError, self).__init__(message + u', caused by ' + repr(cause))
            self.cause = cause
        else:
            super(StructureError, self).__init__(message)


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
        self.forcefield = None

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

        num_chains = self.topology.getNumChains()

        if ends is not None:
            num_ends = len(ends)
            if num_ends != num_chains:
                emsg = 'Number of terminal capping groups ({}) != '.format(num_ends)
                emsg += 'number of chains in the molecule ({})'.format(num_chains)
                raise StructureError(emsg)

            for chain_idx, chain in enumerate(ends, start=1):
                n_cap, c_cap = chain
                if n_cap not in _allowed_n_caps:
                    emsg = 'User-specified N-terminal for chain {}: {}'.format(chain_idx, n_cap)
                    raise StructureError(emsg)

                if c_cap not in _allowed_c_caps:
                    emsg = 'User-specified C-terminal for chain {}: {}'.format(chain_idx, c_cap)
                    raise StructureError(emsg)

        else:
            ends = [('ACE', 'NME') for _ in range(num_chains)]

        # PDBFixer is picky with mmCIF files and OpenMM does not
        # export them properly. For now, we save as PDB and re-read.
        with tempfile.TemporaryFile() as handle:
            app.PDBFile.writeFile(self.topology, self.positions, handle)
            handle.seek(0)  # rewind
            s = pf.PDBFixer(pdbfile=handle)

        # Find missing atoms to exclude from termini addition
        # Avoids adding other atoms when calling `add_termini`
        # Clunky but explicit is better than implicit!
        s.findMissingResidues()
        s.findMissingAtoms()
        missing_residues = set(s.missingResidues.keys())
        missing_atoms = set(s.missingAtoms.keys())

        # Manually build Sequence object(s) and add to structure
        sequences = []
        for chain_idx, chain in enumerate(s.topology.chains()):
            chain_reslist = [r.name for r in chain.residues()]
            n_ter, c_ter = chain_reslist[0], chain_reslist[-1]

            # Add caps if necessary
            n_cap, c_cap = ends[chain_idx]
            if n_cap and n_ter != n_cap:
                chain_reslist.insert(0, n_cap)
            if c_cap and c_ter != c_cap:
                chain_reslist.append(c_cap)

            sequences.append(Sequence(chain.id, chain_reslist))

        s.sequences = sequences
        self.sequences = sequences

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
