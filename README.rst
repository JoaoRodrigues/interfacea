interfacea
=====================================

.. image:: https://zenodo.org/badge/136096537.svg
   :target: https://zenodo.org/badge/latestdoi/136096537

.. image:: https://img.shields.io/travis/JoaoRodrigues/interfacea/master.svg?label=TravisCI
    :target: https://travis-ci.com/JoaoRodrigues/interfacea
    :alt: Travis Continuous Integration

.. image:: https://img.shields.io/appveyor/ci/joaorodrigues/interfacea?label=Appveyor
  :target: https://ci.appveyor.com/project/JoaoRodrigues/interfacea
  :alt: Appveyor Continuous Integration

.. image:: https://api.codeclimate.com/v1/badges/ca054443ee84f96d748a/test_coverage
   :target: https://codeclimate.com/github/JoaoRodrigues/interfacea/test_coverage
   :alt: Test Coverage

.. start-description

interfacea is a Python library to analyze geometric and energetic features of
protein interfaces. The feature list includes:

* Fast and flexible graph-based search of functional groups.

* Identification of several types of atomic interactions
  (e.g. ionic, hydrogen bonds).

* Rebuilding missing (heavy) atoms and *in silico* mutagenesis.

.. end-description

.. start-intro

Getting Started
-----------------

**Installation instructions**

.. code-block:: bash

    git clone https://github.com/joaorodrigues/interfacea.git interfacea
    cd interfacea
    conda env create -f requirements.yml
    python setup.py build && python setup.py install

**Quick Example(s)**

.. code-block:: python

    # Load a molecule and perform energy minimization.

    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')
    mol.prepare(minimize=True)

    print(mol.potential_energy)

.. code-block:: python

    # Find ionic interactions in the structure

    analyzer = ia.InteractionAnalyzer(mol)
    analyzer.get_ionic()

    print(analyzer.itable._table)  # will obviously change in the future.

.. end-intro

Read the full documentation at `interfacea.readthedocs.io <https://interfacea.readthedocs.io/en/latest/>`_

Software Dependencies
---------------------

interfacea requires Python 3.6 or a more recent version, as well as the following packages:

- networkx (https://networkx.github.io)
- pandas (https://http://pandas.pydata.org)
- PDBFixer (https://github.com/pandegroup/pdbfixer)
- OpenMM (http://openmm.org)

These dependencies can be installed by following the instructions on their
respective websites or by typing ``conda env create -f requirements.yml``.
