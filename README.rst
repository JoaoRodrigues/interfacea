interfacea
======================================

.. image:: https://landscape.io/github/JoaoRodrigues/interfacea/master/landscape.svg?style=flat
   :target: https://landscape.io/github/JoaoRodrigues/interfacea/master
   :alt: Landscape.io Code Health

.. image:: https://ci.appveyor.com/api/projects/status/tcss5si0bgtdl3xj?svg=true
   :target: https://ci.appveyor.com/project/JoaoRodrigues/interfacea
   :alt: Appveyor Continuous Integration

.. image:: https://travis-ci.com/JoaoRodrigues/interfacea.svg?branch=master
   :target: https://travis-ci.com/JoaoRodrigues/interfacea
   :alt: Travis Continuous Integration

interfacea (interface energy analysis) is a Python library written to facilitate the analysis of
the energetics of protein interactions. It leverages the high-performance of OpenMM to calculate
pairwise energies between atoms, residues, or chains under popular force fields.

Quick Example
=============

::

    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    mol.add_termini()
    mol.add_missing_atoms()


Software Dependencies
=====================

interfacea runs on **Python 3.x** only and depends on the following packages:
- `OpenMM <http://openmm.org>`
- `PDBFixer <https://github.com/pandegroup/pdbfixer>`

Dependencies can be installed via ``conda`` following instructions on their
websites.