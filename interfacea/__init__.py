#!/usr/bin/env python

"""
Code to calculate pairwise residue energies in macromolecular
structures.
"""

import logging
import os

import simtk.openmm.app as app

from .structure import Structure
from .structure import StructureError


# Setup logger
# This is the parent logger since the library is supposed
# to be loaded from here. Hence, configs set here apply to
# all module-level loggers
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Methods


def read(fpath, ftype=None):
    """Creates a `Structure` instance from a PDB/mmCIF file.

    Args:
        fpath (str): path to file to be parsed.
        ftype (:obj:`str`, optional): file type (PDB or mmCIF). `None` defaults
            to guessing the file type from the extension.

    Returns:
        :obj:`Structure`: new instance of `Structure` class defining a topology
            and positions for the input structure.

    Raises:
        StructureError: generic error class for problems during structure parsing
            or conversion with OpenMM.
    """

    _pdb_formats = {'pdb', 'ent'}
    _cif_formats = {'cif', 'mmcif'}
    _formats = _pdb_formats | _cif_formats
    _formats_str = ','.join(_formats)

    logger.debug('Reading file: {}'.format(fpath))

    fullpath = os.path.abspath(fpath)
    if not os.path.isfile(fullpath):
        emsg = 'File could not be read: {}'.format(fullpath)
        raise StructureError(emsg)

    fname, fext = os.path.splitext(fullpath)
    ftype = ftype if ftype is not None else fext[1:]
    logger.debug('Assigning file type: {}'.format(ftype))

    if ftype in _pdb_formats:
        try:
            struct = app.PDBFile(fullpath)
            logger.debug('Successfully parsed file as PDB')
        except Exception as e:
            emsg = 'Failed parsing file {} as \'PDB\' format'.format(fullpath)
            raise StructureError(emsg, e)

    elif ftype in _cif_formats:
        try:
            struct = app.PDBxFile(fullpath)
            logger.debug('Successfully parsed file as mmCIF')
        except Exception as e:
            emsg = 'Failed parsing file {} as \'mmCIF\' format'.format(fullpath)
            raise StructureError(emsg, e)
    else:
        emsg = '\'{}\' is not one of the supported types: {}'.format(ftype, _formats_str)
        raise StructureError(emsg)

    return Structure(fullpath, struct)
