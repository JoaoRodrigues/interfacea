
interfacea
=====================================

.. image:: https://img.shields.io/travis/JoaoRodrigues/interfacea/master.svg?label=TravisCI
    :target: https://travis-ci.com/JoaoRodrigues/interfacea
    :alt: Travis Continuous Integration

.. image:: https://img.shields.io/appveyor/ci/joaorodrigues/interfacea?label=Appveyor
  :target: https://ci.appveyor.com/project/JoaoRodrigues/interfacea
  :alt: Appveyor Continuous Integration

interfacea is a Python library to analyze geometric and energetic features of
protein interfaces. The feature list includes:

* Fast and flexible graph-based search of functional groups.

* Identification of several types of atomic interactions
  (e.g. ionic, hydrogen bonds).

* Rebuilding missing (heavy) atoms and *in silico* mutagenesis.

For the Impatient
-----------------

Read the full documentation at `interfacea.readthedocs.io <https://interfacea.readthedocs.io/en/latest/>`_.

**Installation instructions**

.. code:: bash

    git clone https://github.com/joaorodrigues/interfacea.git interfacea
    cd interfacea
    conda env create -f requirements.yml  # create environment with dependencies
    python setup.py build && python setup.py install  # install interfacea

**Quick Example(s)**

* Loading a molecule and performing energy minimization.

    .. code-block:: python

        import interfacea as ia

        mol = ia.read('tests/data/mini.pdb')
        mol.prepare(minimize=True)

        print(mol.potential_energy)

* Find ionic interactions in the structure

    .. code-block:: python

        import interfacea as ia

        mol = ia.read('tests/data/mini.pdb')

        analyzer = ia.InteractionAnalyzer(mol)
        analyzer.get_ionic()

        print(analyzer.itable._table)  # will obviously change in the future.


Software Dependencies
---------------------

interfacea requires Python 3.6 or a more recent version, as well as the following packages:

- networkx (https://networkx.github.io)
- pandas (https://http://pandas.pydata.org)
- PDBFixer (https://github.com/pandegroup/pdbfixer)
- OpenMM (http://openmm.org)

These dependencies can be installed by following the instructions on their
respective websites or by typing ``conda env create -f requirements.yml``.
