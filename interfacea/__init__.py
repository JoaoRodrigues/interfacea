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
import pathlib
import random

# Setup logger
# This is the parent logger since the library is supposed
# to be loaded from here. Hence, configs set here apply to
# all module-level loggers
logging.getLogger(__name__).addHandler(logging.NullHandler())
logging.getLogger(__name__).setLevel(logging.CRITICAL)

# Global Constants
RANDOM_SEED = random.randint(0, 1000)  # user can set it manually later


# Methods
def set_log_level(level='verbose'):
    """Enables logging to a certain level.

    Useful for interactive/debugging applications.

    Args:
        level (str): verbosity/type of logging. Choose from:
            - 'silent': disables logging.
            - 'minimal': only warnings and other critical messages
            - 'verbose': informative/descriptive messages (default).
            - 'debug': very verbose internal/diagnostic messages.
    """

    # Logging levels
    _lc = {
           'minimal': logging.WARNING,
           'verbose': logging.INFO,
           'debug': logging.DEBUG
          }

    if level not in _lc and level != 'silent':
        emsg = f"Unknown or unsupported log level: {level}"
        raise ValueError(emsg)

    # Treat 'silent' differently
    if level == 'silent':
        root_logger = logging.getLogger()
        root_logger.handlers = []  # clear handler list
        root_logger.setLevel(logging.CRITICAL)
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

    log_level = _lc.get(level)
    root_logger.setLevel(log_level)
    logging.warn(f"Logging activated and set to '{level}'")  # always show


# Randomness
def set_random_seed(seed=917):
    """Sets a defined seed for reproducible operations across the library.

    This does not ensure *complete reproducibility*. Some methods in OpenMM,
    for example, are not deterministic across different hardware configurations
    even with the same random seed.
    """

    global RANDOM_SEED

    if isinstance(seed, int) and seed > 0:
        RANDOM_SEED = seed
    else:
        emsg = f"Random seed must be a positive integer: {seed}"
        raise TypeError(emsg)


# IO
def read(fpath, ftype=None, **kwargs):
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

    logging.info(f"Reading file: {fpath}")

    fpath = pathlib.Path(fpath)

    if not fpath.exists():
        emsg = f"File unreadable/does not exist: {fpath}"
        raise StructureError(emsg)

    fext = fpath.suffix[1:]
    if ftype is None:
        ftype = fext

    logging.debug(f'File type: {ftype}')

    if ftype in _pdb_formats:
        try:
            struct = app.PDBFile(str(fpath.resolve()))
        except Exception as e:
            emsg = f"Could not parse file as PDB: {fpath}"
            raise StructureError(emsg) from e

    elif ftype in _cif_formats:
        try:
            struct = app.PDBxFile(str(fpath.resolve()))
        except Exception as e:
            emsg = f"Could not parse file as mmCIF: {fpath}"
            raise StructureError(emsg) from e

    else:
        emsg = f"Unsupported file type: {ftype}"
        raise StructureError(emsg)

    return Structure(fpath.name, struct, **kwargs)
