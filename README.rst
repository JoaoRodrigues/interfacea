interfacea
======================================

.. image:: https://landscape.io/github/JoaoRodrigues/interfacea/master/landscape.svg?style=flat
   :target: https://landscape.io/github/JoaoRodrigues/interfacea/master
   :alt: Landscape.io Code Health

.. image:: https://ci.appveyor.com/api/projects/status/tcss5si0bgtdl3xj?svg=true
   :target: https://ci.appveyor.com/project/JoaoRodrigues/interfacea
   :alt: Appveyor Continuous Integration


interfacea (interface energy analysis) is a Python library written to facilitate the analysis of
the energetics of protein interactions. It leverages the high-performance of OpenMM to calculate
pairwise energies between atoms, residues, or chains under popular force fields.

Quick Example
=============

To use interfacea you need OpenMM (see http://openmm.org) and PDBFixer (see https://github.com/pandegroup/pdbfixer/) 
installed ::

    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    mol.add_termini()
    mol.add_missing_atoms()


