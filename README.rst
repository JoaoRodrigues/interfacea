interfacea
=====================================

.. image:: https://zenodo.org/badge/136096537.svg
   :target: https://zenodo.org/badge/latestdoi/136096537

.. image:: https://api.codeclimate.com/v1/badges/ca054443ee84f96d748a/test_coverage
   :target: https://codeclimate.com/github/JoaoRodrigues/interfacea/test_coverage
   :alt: Test Coverage

.. image:: https://dev.azure.com/jpglmrodrigues/interfacea-ci/_apis/build/status/JoaoRodrigues.interfacea?branchName=reorg&label=Build
   :target: https://dev.azure.com/jpglmrodrigues/interfacea-ci/_build/latest?definitionId=1&branchName=reorg
   :alt: Azure Build Check

.. start-description

interfacea is a Python library to analyze geometric and energetic features of
protein interfaces. It implements a graph-based functional group
representation that allows for fast detection of several types of chemical
interactions. In addition, interfacea is tightly coupled with
`OpenMM <http://openmm.org/>`_ to provide simple modeling and refinement
operations on atomic structures under a variety of force fields.

We welcome contributions and suggestions to our project - please read our
contributing guidelines.

.. end-description

.. start-intro

Getting Started
-----------------

**Installing from Source**

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/joaorodrigues/interfacea.git interfacea

    cd interfacea

    # Setup a dedicated conda environment to install interfacea
    conda create --yes --quiet --name interfacea python=3.6
    conda env update --name interfacea --file environment.yml

    # Activate the environment and install the library
    conda activate interfacea
    python setup.py build
    python setup.py install

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
