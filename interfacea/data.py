#!/usr/bin/env python

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
