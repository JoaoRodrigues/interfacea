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
import random

__all__ = ['set_random_seed', 'set_log_level']

logging.getLogger(__name__).addHandler(logging.NullHandler())

# Constants
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
