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

from .interactions import InteractionAnalyzer

# Setup logger
# This is the parent logger since the library is supposed
# to be loaded from here. Hence, configs set here apply to
# all module-level loggers
logging.getLogger(__name__).addHandler(logging.NullHandler())


# Methods
def set_log_level(level='none'):
    """Enables logging to a certain level.

    Useful for interactive/debugging applications.

    Args:
        level (str): verbosity/type of logging. Can be either
            'none', 'minimal', or 'verbose'.
    """

    if level == 'none':
        pass
    elif level == 'minimal':
        # Add timestamp
        logging.basicConfig(format='[%(asctime)s] %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.INFO)
    elif level == 'verbose':
        logging.basicConfig(format='[%(asctime)s] %(module)s :: %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
    else:
        raise ValueError('Logging level must be: \'none\', \'minimal\', or \'verbose\'')

    logging.info('Logging enable and set to \'{}\''.format(level))


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

    logging.debug('Reading file: {}'.format(fpath))

    fullpath = os.path.abspath(fpath)
    if not os.path.isfile(fullpath):
        emsg = 'File could not be read: {}'.format(fullpath)
        raise StructureError(emsg)

    fname, fext = os.path.splitext(fullpath)
    ftype = ftype if ftype is not None else fext[1:]
    logging.debug('Assigning file type: {}'.format(ftype))

    if ftype in _pdb_formats:
        try:
            logging.debug('Detected file as PDB')
            struct = app.PDBFile(fullpath)
        except Exception as e:
            emsg = 'Failed parsing file {} as \'PDB\' format'.format(fullpath)
            raise StructureError(emsg) from e

    elif ftype in _cif_formats:
        try:
            logging.debug('Detected file as mmCIF')
            struct = app.PDBxFile(fullpath)
        except Exception as e:
            emsg = 'Failed parsing file {} as \'mmCIF\' format'.format(fullpath)
            raise StructureError(emsg) from e
    else:
        emsg = '\'{}\' is not one of the supported types: {}'.format(ftype, _formats_str)
        raise StructureError(emsg)

    logging.info('Successfully parsed file into Structure: {}'.format(fullpath))
    return Structure(fullpath, struct)
