#!/usr/bin/env python

"""
Unit tests for read() function
"""

import os
import sys
import unittest

if __name__ == '__main__':

    mpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    sys.path.insert(0, mpath)  # module path so we load dist files

    loader = unittest.defaultTestLoader

    tpath = os.path.join(mpath, 'tests')
    suite = loader.discover(tpath)

    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
