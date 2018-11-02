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
Analysis of Biomolecular Interfaces.
"""

from __future__ import print_function

import logging

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())

# Protein amino acid names (included substitutions)
# Taken from pdbfixer.py
protein_aa = {'2AS', '3AH', '5HP', 'ACL', 'AGM', 'AIB', 'ALA', 'ALM', 'ALO',
              'ALY', 'ARG', 'ARM', 'ASA', 'ASB', 'ASK', 'ASL', 'ASN', 'ASP',
              'ASQ', 'AYA', 'BCS', 'BHD', 'BMT', 'BNN', 'BUC', 'BUG', 'C5C',
              'C6C', 'CAS', 'CCS', 'CEA', 'CGU', 'CHG', 'CLE', 'CME', 'CSD',
              'CSO', 'CSP', 'CSS', 'CSW', 'CSX', 'CXM', 'CY1', 'CY3', 'CYG',
              'CYM', 'CYQ', 'CYS', 'DAH', 'DAL', 'DAR', 'DAS', 'DCY', 'DGL',
              'DGN', 'DHA', 'DHI', 'DIL', 'DIV', 'DLE', 'DLY', 'DNP', 'DPN',
              'DPR', 'DSN', 'DSP', 'DTH', 'DTR', 'DTY', 'DVA', 'EFC', 'FLA',
              'FME', 'GGL', 'GL3', 'GLN', 'GLU', 'GLY', 'GLZ', 'GMA', 'GSC',
              'HAC', 'HAR', 'HIC', 'HIP', 'HIS', 'HMR', 'HPQ', 'HTR', 'HYP',
              'IAS', 'IIL', 'ILE', 'IYR', 'KCX', 'LEU', 'LLP', 'LLY', 'LTR',
              'LYM', 'LYS', 'LYZ', 'MAA', 'MEN', 'MET', 'MHS', 'MIS', 'MLE',
              'MPQ', 'MSA', 'MSE', 'MVA', 'NEM', 'NEP', 'NLE', 'NLN', 'NLP',
              'NMC', 'OAS', 'OCS', 'OMT', 'PAQ', 'PCA', 'PEC', 'PHE', 'PHI',
              'PHL', 'PR3', 'PRO', 'PRR', 'PTR', 'PYX', 'SAC', 'SAR', 'SCH',
              'SCS', 'SCY', 'SEL', 'SEP', 'SER', 'SET', 'SHC', 'SHR', 'SMC',
              'SOC', 'STY', 'SVA', 'THR', 'TIH', 'TPL', 'TPO', 'TPQ', 'TRG',
              'TRO', 'TRP', 'TYB', 'TYI', 'TYQ', 'TYR', 'TYS', 'TYY', 'VAL'}

# van der Waals radii (in Angstrom) from:
#    A. Bondi, J. Phys. Chem., 68, 441 - 452, 1964,
# radius of H taken from:
#    R.S. Rowland & R. Taylor, J.Phys.Chem., 100, 7384 - 7391, 1996.
vdw_radii = {
    1: 1.20,
    5: 2.00,  # default value
    6: 1.70,
    7: 1.55,
    8: 1.52,
    9: 1.47,
    14: 2.10,
    15: 1.80,
    16: 1.80,
    17: 2.27,
    19: 1.76,
    20: 1.37,
    26: 2.00,  # default
    27: 2.00,  # default
    28: 1.63,
    29: 1.40,
    30: 1.39,
    33: 1.85,
    34: 1.90,
    35: 1.85,
    48: 1.58,
    53: 1.98,
    79: 1.66,
    80: 1.55,
}

# Covalent atom radii (in Angstrom) used to infer connectivity
# Values taken from CCDC and Roger Sayle's website:
# http://www.daylight.com/meetings/mug01/Sayle/m4xbondage.html
covalent_radii = {
    1: 0.23,
    5: 0.83,
    6: 0.68,
    7: 0.68,
    8: 0.68,
    9: 0.64,
    14: 1.20,
    15: 1.05,
    16: 1.02,
    17: 0.99,
    33: 1.21,
    34: 1.22,
    35: 1.21,
    52: 1.47,
    53: 1.40
}
