.. _contrib-label:

Contributing Guidelines
=======================

First off, thanks for being interested in contributing!

This document is a set of guidelines and instructions for contributing
to interfacea and its development. These are guidelines, *not rules*, so
use your best judgement and feel free to propose changes to this document
in a pull request.

Table of Contents
-----------------


1. `Our Code of Conduct`_

2. `What Should I Know Before Getting Started?`_

3. `How Can I Contribute?`_
    - Reporting Bugs
    - Suggesting New Features
    - Contributing Code
    - `Contributing Documentation`_

4. `Style Guides`_
    - Coding Style Guide
    - Documentation Style Guide
    - Git(hub) Etiquette

5. `Continuous Integration`_

Our Code of Conduct
-------------------

This project and everyone participating in it is governed by our
:ref:`Code of Conduct <coc-label>`. By participating, we
expect you to abide by and uphold this code. Please report any
unacceptable behavior to joaor@stanford.edu.

What Should I Know Before Getting Started?
------------------------------------------

For now, nothing really...

How Can I Contribute?
---------------------

Reporting Bugs
~~~~~~~~~~~~~~

Suggesting New Features
~~~~~~~~~~~~~~~~~~~~~~~

Contributing Code
~~~~~~~~~~~~~~~~~

If you are reading this section, thank you for donating your time
to improve the code in interfacea. The following guidelines should
help you get up to speed with our development practices.

**Feature Branches**
We make heavy use of git and github (obviously) to develop interfacea.
The :code:`master` branch holds the latest version of the code. In order
to avoid conflicts during development, we do not accept direct commits to
the :code:`master` branch. Instead, we make use of *feature branches* and
pull requests. In short, assuming you have cloned/fetched the latest
version of the master branch, the recommended pipeline for development is:

.. code-block:: bash

    # Make a new branch::
    git checkout -B my-new-feature

    # Do your magic.
    ...

    # Push the changes to your fork
    git push <your fork> my-new-feature

    # On GitHub, issue a pull request of your fork vs our master branch.

**Development Environment**
Every developer has their own style and preferences. Our aim in writing
this section is not to *impose* on your habits, but to raise awareness
to ours, so that you know what we usually expect in code contributions.

The :code:`devtools` folder contains our recommended development environment,
specified in a conda YAML file. We use `flake8 <http://flake8.pycqa.org>`_ to
assert a coherent style in our Python code - see our `Coding Style Guide`_ for
more details on this. Then, we use `pytest <https://docs.pytest.org/>`_ to run
unit tests and make sure the code runs as expected, and `coverage.py <https://coverage.readthedocs.io/>`_
to quantify how much of the codebase is being tested.

When you issue a pull request to the master branch, our continuous integration
services will run checkers for all three: Python PEP8 compliance via flake8
and unit tests and their coverage via pytest and coverage.py. To make your
life easier, and ours too, we ask you to run (some of these) checks locally on
your development machine before issuing pull requests. To this end, we provide
pre-commit hooks for git to run flake8 via `pre-commit <https://pre-commit.com/>`_
that you are welcome to install on your machine. These obliviate the need to
manually run flake8 every time you want to push, as well as the usual commits
to fix style issues. We **do not** include pytest on these pre-commit hooks
on purpose, as we believe it might slow down development. We leave this up to
you: you can run pytest locally before pushing, or just wait for our CI
servers to run them for you.

To setup your development environment as per our recommendations, issue the
following commands:

.. code-block:: bash

    # Create the development environment
    conda env create -f devtools/dev-environment.yml

    # Run the pre-commit hooks once to setup
    pre-commit run --all-files

    # (Optional)
    # Install the pre-commit hooks so that they run
    # every time you issue git commit
    pre-commit install

    # Develop away...

    # ADD TESTS HERE
    # ....

Contributing Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Thank you for taking the time to contribute to our documentation! We believe
that proper documentation is essential for the widespread adoption and usage
of software tools.

The first step to contribute to our project is to fork our repository,
clone it locally, and create a new branch named after your contribution.
For example, if you want to fix a typo on the documentation index page,
name your branch something like :code:`fix-index-typo` and make your
changes there. Then, once you are done, push it to your github account
and make a pull request to our repository (see `Git(hub) Etiquette`_).

We use `Sphinx <http://www.sphinx-doc.org/en/master/>`_ to
build and format our online documentation, which is then hosted at
`our ReadTheDocs webpage <https://interfacea.readthedocs.io>`_. Every module,
class, and function in interfacea should have docstrings, following the
`Google docstring style guide <http://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings>`_.
These docstrings are then collected and automatically formatted by Sphinx
to build our API reference pages. Other pages, like this one, are manually
write in Restructured Text format (see a syntax guide
`here <https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_).

If you were to add a new module and wanted it to show in the API reference
pages, you would start by creating a document in the :code:`docs/api-ref`
folder. This file should have a :code:`.. currentmodule::` directive
pointing to your new code and an :code:`.. autosummary::` directive
followed by the classes or functions you wish to feature in the
documentation. Have a look at :code:`docs/api-ref/interfacea.core.rst`
for a starting template. Feel free to add more text to the page(s) you
create, to contextualize your code and to explain its logic. If
appropriate, you should also add an *Examples* section showcasing your
code, or just a blurb in the :code:`docs/examples.rst` page.

If you find an error, omission, or think something could be phrased more
clearly, please start by reporting it in our issue tracker on GitHub. Also,
read the following paragraphs as they try to explain (roughly) how the
documentation is built and which files of the source code influence that
process.

* The :code:`docs/` folder contains the documents that will be rendered by
  Sphinx and ReadTheDocs. The makefiles are whatever Sphinx left when we
  ran it the first time and can be used to preview how the documentation
  will look once online.

* You should use the :code:`docs/docs-environment.yml` file with conda-env
  (:code:`conda env create -f docs-environment.yml`) to create an environment
  named *interfacea-docs* with all the dependencies necessary to build and
  test the documentation. Then, in the :code:`docs/` folder, issue
  :code:`make html` to build the documentation. You can preview it by
  opening the newly-created :code:`docs/_build/html/index.html` file in
  your web browser.

* The :code:`docs/conf.py` file contains the configuration for the Sphinx
  build. This file has settings such as the theme we chose to render the
  pages, which plugins we make use of, and other small but important details.

* The last piece of the documentation puzzle is the configuration file for
  ReadTheDocs: :code:`.readthedocs.yml`. Unlike the rest of the
  documentation, this file lives at the root of the source code repository
  and is used by the ReadTheDocs server to start our online build.

Style Guides
------------

Coding Style Guide
~~~~~~~~~~~~~~~~~~

The code in interfacea *must* be written in Python 3.x and follow the PEP8.
Our flake8 linter has one ignore rule, about lambda functions, and that's
pretty much it. We use `snake_case` instead of `CamelCase` for variables
and methods, and `CapWords` for class names. We specifically *do not*
use `black <https://github.com/psf/black>`_ because we find it too restrictive.

Documentation Style Guide
~~~~~~~~~~~~~~~~~~~~~~~~~

* Use `Restructured Text format <https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
* Use the `Google conventions for docstrings <http://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings>`_

Git(hub) Etiquette
~~~~~~~~~~~~~~~~~~

Please try to follow these guidelines when writing commit messages:
* Write a comprehensive message: no `Fixed bug` messages please!
* Limit the first line to 72 characters.
* Reference issues and pull requests after the first line.

# TODO
# Add what labels we use and what they are for.
# Both issue and PR labels.

Continuous Integration
----------------------

We aim to have all the building, testing, and documentation infrastructure
automated through continuous integration services/servers. In this way, we
hope to make it easier for new and current developers alike to contribute
and maintain the library. In accordance with the previous sections in this
document, we have several checks that run automatically with each pull
request that is made against the master branch. The following paragraphs
explain, roughly, how this process happens and the choices we made in terms
to host/carry out the continuous integration services.

The bulk of the services run on `Azure Pipelines <https://dev.azure.com/jpglmrodrigues/interfacea-ci>`_.
We have pipelines setup to 1) ensure if the library builds correctly,
2) run the testing framework (and quantify how much of the code they
cover), and 3) verify if the documentation builds correctly. Each of this
actions is setup in a *stage* in Azure DevOps terminology.

These pipelines are controlled by the :code:`azure-pipelines.yml` file, as
well as its dependencies in the :code:`ci/` folder. By default, we run
builds and the testing framework against several versions of Python - using
tox - and all three major OSes (Windows, Mac OSX, and Linux).

As with the rest of the library, we welcome pull requests that improve or
correct our continuous integration practices!
