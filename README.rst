
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


Installation
============

```bash
git clone https://github.com/joaorodrigues/interfacea.git interfacea
cd interfacea
conda env create -f requirements.yml  # create environment with dependencies
python setup.py build && python setup.py install  # install interfacea
```


Quick Example(s)
=============

1. Loading a molecule and performing energy minimization.

    ```
    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    mol.prepare(minimize=True)

    print(mol.potential_energy)
    ```

2. Find ionic interactions in the structure

    ```
    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    
    analyzer = ia.InteractionAnalyzer(mol)
    analyzer.get_ionic()

    print(analyzer.itable._table)  # will obviously change in the future.
    ```

3. Find specific functional groups in the structure

    ```
    import interfacea as ia
    import interfacea.functional_groups as fgs

    mol = ia.read('tests/data/mini.pdb')
    
    # Find phosphate groups in molecule

    phosphate = fgs.Phosphate()
    matches = phosphate.search(mol)
    for res, groups in matches.items:
        n_groups = len(groups)
        if n_groups:
            print('Residue {}{} contains {} phosphate(s)'.format(res.name, res.id, n_groups))
    ```

Software Dependencies
=====================

interfacea runs on **Python 3.6+** only and requires the following packages:

- networkx (https://networkx.github.io)
- pandas (https://http://pandas.pydata.org)
- PDBFixer (https://github.com/pandegroup/pdbfixer)
- OpenMM (http://openmm.org)

Dependencies should be installed via ``conda`` using
``conda env create -f requirements.yml`` or by following instructions on their
websites.
