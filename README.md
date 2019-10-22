
[![Travis Continuous Integration](https://travis-ci.com/JoaoRodrigues/interfacea.svg?branch=master)](https://travis-ci.com/JoaoRodrigues/interfacea)
[![Appveyor Continuous Integration](https://ci.appveyor.com/api/projects/status/tcss5si0bgtdl3xj?svg=true)](https://ci.appveyor.com/project/JoaoRodrigues/interfacea)
[![DOI](https://zenodo.org/badge/136096537.svg)](https://zenodo.org/badge/latestdoi/136096537)



interfacea
======================================

interfacea is a Python library written to facilitate the analysis of interfaces
of a wide range of biomolecular complexes. The feature list includes:
  * Fast and flexible graph-based search of functional groups.
  * Identification of several types of atomic interactions (e.g. ionic, hydrogen bonds).
  * Rebuilding missing (heavy) atoms and _in silico_ mutagenesis.
  * Calculation and decomposition of energetics using OpenMM.
  * ...


Installation
------------

```bash
git clone https://github.com/joaorodrigues/interfacea.git interfacea
cd interfacea
conda env create -f requirements.yml  # create environment with dependencies
python setup.py build && python setup.py install  # install interfacea
conda activate interfacea-install
# Ta-Daa!
```


Quick Example(s)
----------------

1. Loading a molecule and performing energy minimization.

    ```python
    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    mol.prepare(minimize=True)

    print(mol.potential_energy)
    ```

2. Find ionic interactions in the structure

    ```python
    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    
    analyzer = ia.InteractionAnalyzer(mol)
    analyzer.get_ionic()

    print(analyzer.itable._table)  # will obviously change in the future.
    ```

3. Find specific functional groups in the structure

    ```python
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
---------------------

interfacea runs on **Python 3.6+** only and requires the following packages:

- networkx (https://networkx.github.io)
- pandas (https://http://pandas.pydata.org)
- PDBFixer (https://github.com/pandegroup/pdbfixer)
- OpenMM (http://openmm.org)

Dependencies should be installed via ``conda`` using
``conda env create -f requirements.yml`` or by following instructions on their
websites.
