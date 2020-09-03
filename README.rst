interfacea
=====================================

.. .. image::
..    :target:
..    :alt: Zenodo DOI

.. .. image::
..    :target:
..    :alt: Test Suite

.. .. image::
..    :target:
..    :alt: Test Coverage

.. start-description

interfacea is a Python library to analyze protein interfaces. You can use it to find
hydrogen bonds and other chemical interactions in your 3D models, as well as perform
(simple) *in silico* mutagenesis.

Our goal is to provide an open-source and extensible framework that enables you both
to interactively explore and analyze your protein structures, as a user, or to use
as a stepping stone to develop your own tools.

We gladly welcome contributions and suggestions to our code and documentation.

.. end-description

.. start-intro

Getting Started
-----------------

You can find the full documentation at
`interfacea.readthedocs.io <https://interfacea.readthedocs.io/en/latest/>`_

**Installing from Source**

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/joaorodrigues/interfacea.git interfacea

    cd interfacea

    # Setup a dedicated conda environment to install interfacea
    conda env create -f environment.yml
    conda activate interfacea

    # Install the library
    python setup.py install

**Quick Example(s)**

.. code-block:: python

    import interfacea as ia

.. end-intro

Software Dependencies
---------------------

interfacea requires a modern version of Python (3.8+) and depends on the following
packages:

- `numpy <https://numpy.org/>`_
- `networkx <https://networkx.github.io>`_

Certain submodules might require additional dependencies. If you are not picky about
disk space (we won't use much anyway), you can simply use our ```environment.yml```
file as shown in the Installing from Source section above.
