#!/usr/bin/env python

"""
Setup script for packaging/installing interfacea.
"""

from distutils.core import setup

setup(
    name='interfacea',
    version='0.0-dev',
    packages=['interfacea', ],
    data_files=[("", ["LICENSE"])],
    license='Apache License',
    long_description="Package to analyze the energetics of biomolecular interfaces.",
)
