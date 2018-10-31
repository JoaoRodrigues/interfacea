#!/usr/bin/env python

"""
Analysis of Biomolecular Interfaces.
"""

from __future__ import print_function

import logging

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())

# van der Waals radii from:
#    A. Bondi, J. Phys. Chem., 68, 441 - 452, 1964,
# radius of H taken from:
#    R.S. Rowland & R. Taylor, J.Phys.Chem., 100, 7384 - 7391, 1996.
vdw_radii = {
    1: 1.20,
    6: 1.70,
    7: 1.55,
    8: 1.52,
    9: 1.47,
    15: 1.80,
    16: 1.80,
}
