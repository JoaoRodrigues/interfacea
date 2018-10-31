#!/usr/bin/env python

"""
Code to calculate pairwise residue energies in macromolecular structures.
"""

import logging

__all__ = ['set_random_seed', 'set_log_level']

logging.getLogger(__name__).addHandler(logging.NullHandler())

# Constants
RANDOM_SEED = None


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
