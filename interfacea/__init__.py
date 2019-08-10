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
Code to calculate pairwise residue energies in macromolecular structures.
"""

import logging
import os
import random

# Setup logger
# This is the parent logger since the library is supposed
# to be loaded from here. Hence, configs set here apply to
# all module-level loggers
logging.getLogger(__name__).addHandler(logging.NullHandler())

# Global Constants
RANDOM_SEED = random.randint(0, 1000)  # user can set it manually later


# Methods
def set_log_level(level='minimal'):
    """Enables logging to a certain level.

    Useful for interactive/debugging applications.

    Args:
        level (str): verbosity/type of logging. Can be either
            'none', 'minimal', or 'verbose'. Default is 'minimal'.
    """

    if level == 'none':
        root_logger = logging.getLogger()
        root_logger.handlers = []  # clear handler list
        root_logger.setLevel(logging.WARNING)
        return

    handler = logging.StreamHandler()
    formatter = logging.Formatter(fmt='[%(asctime)s] %(message)s',
                                  datefmt='%H:%M:%S')
    handler.setFormatter(formatter)

    # We override the root logger here, assuming this function is only called
    # interactively ...
    root_logger = logging.getLogger()
    root_logger.handlers = []  # clear handler list
    root_logger.addHandler(handler)

    if level == 'minimal':
        root_logger.setLevel(logging.INFO)
    elif level == 'verbose':
        root_logger.setLevel(logging.DEBUG)
    else:
        raise ValueError('Logging level must be: \'none\', \'minimal\', or \'verbose\'')

    logging.info('Logging enabled and set to \'{}\''.format(level))


# Randomness
def set_random_seed(seed=917):
    """Sets a defined seed for reproducible operations across the library.

    This does not ensure *complete reproducibility*. Some methods in OpenMM, for
    example, are not deterministic across different hardware configurations even
    with the same random seed.
    """

    global RANDOM_SEED

    if isinstance(seed, int):
        RANDOM_SEED = seed
    else:
        emsg = 'Invalid random seed: {} - Must be a positive integer.'
        raise TypeError(emsg.format(seed))


# IO
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

    import simtk.openmm.app as app
    from interfacea.core.Structure import Structure
    from interfacea.core.Structure import StructureError

    _pdb_formats = {'pdb', 'ent'}
    _cif_formats = {'cif', 'mmcif'}
    _formats = _pdb_formats | _cif_formats
    _formats_str = ','.join(_formats)

    logging.info('Reading file: {}'.format(fpath))

    fullpath = os.path.abspath(fpath)
    if not os.path.isfile(fullpath):
        emsg = 'File could not be read: {}'.format(fullpath)
        raise StructureError(emsg)

    fname, fext = os.path.splitext(fullpath)
    ftype = ftype if ftype is not None else fext[1:]
    logging.debug('Assigned file type: {}'.format(ftype))

    if ftype in _pdb_formats:
        try:
            struct = app.PDBFile(fullpath)
        except Exception as e:
            emsg = 'Failed parsing file {} as \'PDB\' format'.format(fullpath)
            raise StructureError(emsg) from e

    elif ftype in _cif_formats:
        try:
            struct = app.PDBxFile(fullpath)
        except Exception as e:
            emsg = 'Failed parsing file {} as \'mmCIF\' format'.format(fullpath)
            raise StructureError(emsg) from e
    else:
        emsg = '\'{}\' is not one of the supported types: {}'.format(ftype, _formats_str)
        raise StructureError(emsg)

    # logging.debug('File parsed successfully using OpenMM reader.')
    return Structure(os.path.basename(fname), struct)
