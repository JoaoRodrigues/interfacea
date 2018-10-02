
.. image:: https://landscape.io/github/JoaoRodrigues/interfacea/master/landscape.svg?style=flat
   :target: https://landscape.io/github/JoaoRodrigues/interfacea/master
   :alt: Landscape.io Code Health

.. image:: https://ci.appveyor.com/api/projects/status/tcss5si0bgtdl3xj?svg=true
   :target: https://ci.appveyor.com/project/JoaoRodrigues/interfacea
   :alt: Appveyor Continuous Integration

.. image:: https://travis-ci.com/JoaoRodrigues/interfacea.svg?branch=master
   :target: https://travis-ci.com/JoaoRodrigues/interfacea
   :alt: Travis Continuous Integration


interfacea
======================================

interfacea (interface energy analysis) is a Python library written to facilitate the analysis of
the energetics of protein interactions. It leverages the high-performance of OpenMM to calculate
pairwise energies between atoms, residues, or chains under popular force fields.

Quick Example(s)
=============

::

    import interfacea as ia

    # Load molecule
    mol = ia.read('tests/data/mini.pdb')
    
    # Add missing atoms and cap termini with ACE/NME groups
    mol.add_termini()
    mol.add_missing_atoms()
    
    # Parameterize the molecule using the Amber14 FF in OpenMM
    # Add protons using the definition of the forcefield.
    mol.parameterize(forcefield='amber14-all.xml')
    mol.protonate(ph=7.0)
    
    # Analyze and categorize residue interactions based on geometric criteria
    analyzer = ia.InteractionAnalyzer(mol)
    analyzer.find_salt_bridges()

    # Find phosphate groups in molecule
    import ia.functional_groups as fgs
    phosphate = fgs.Phosphate()
    for residue in mol.topology.residues():
        if phosphate.match(residue):
            print('Residue {} contains phosphate(s)')


Software Dependencies
=====================

interfacea runs on **Python 3.x** only and depends on the following packages:

- OpenMM (<http://openmm.org>
- PDBFixer (<https://github.com/pandegroup/pdbfixer>)
- networkx (<https://networkx.github.io>)

Dependencies can be installed via ``conda`` following instructions on their
websites.
