interfacea
=====================================

.. image:: https://zenodo.org/badge/136096537.svg
   :target: https://zenodo.org/badge/latestdoi/136096537
   :alt: Zenodo DOI

.. image:: https://dev.azure.com/jpglmrodrigues/interfacea-ci/_apis/build/status/JoaoRodrigues.interfacea?branchName=reorg&label=Build
   :target: https://dev.azure.com/jpglmrodrigues/interfacea-ci/_build/latest?definitionId=1&branchName=reorg
   :alt: Azure Build Check

.. image:: https://codecov.io/gh/JoaoRodrigues/interfacea/branch/refactor_v1/graph/badge.svg
   :target: https://codecov.io/gh/JoaoRodrigues/interfacea
   :alt: Test Coverage

.. image:: https://img.shields.io/lgtm/alerts/g/JoaoRodrigues/interfacea.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/JoaoRodrigues/interfacea/alerts/
   :alt: LGTM Maintainability Score

.. start-description

**Disclaimer**
Our code is undergoing a drastic remodelling - we do not guarantee it works! Look at the master branch for a working (outdated) version.

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
    conda create --yes --quiet --name interfacea python=3.7
    conda env update --name interfacea --file environment.yml

    # Activate the environment and install the library
    conda activate interfacea
    python setup.py build
    python setup.py install

**Quick Example(s)**

.. code-block:: python

    # Load a molecule.

    import interfacea as ia

    mol = ia.read('tests/data/mini.pdb')

.. end-intro

Read the full documentation at `interfacea.readthedocs.io <https://interfacea.readthedocs.io/en/latest/>`_

Software Dependencies
---------------------

interfacea requires Python 3.7 and depends on the following packages:

- numpy (https://numpy.org/)
- networkx (https://networkx.github.io)

You can install these dependencies manually or let us take care of it by
doing ``conda env create -f requirements.yml``.
