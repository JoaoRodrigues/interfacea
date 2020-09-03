.. _contrib-label:

Contributing Guidelines
=======================

First of all, thank you for being interested in contributing! We really appreciate
your willingness to donate your time and skills to improving interfacea.

This document puts forth a set of guidelines for contributing to interfacea.
Importantly, these are *guidelines* and not *rules*. Use your best judgement and do
feel free to propose changed to this very document in a pull request!

Also, make sure you read our :ref:`Code of Conduct <coc-label>`. By participating in
the development of interfacea, we expect you to abide by and enforce this code.
Please report any unacceptable behavior to the maintainers (e.g. @joaorodrigues).

Table of Contents
-----------------

1. `How Can I Contribute?`_
    - Reporting Bugs
    - Suggesting New Features
    - Contributing Code

2. `Style Guides`_
    - Coding Style Guide
    - Git(hub) Etiquette


How Can I Contribute?
---------------------

Reporting Bugs
~~~~~~~~~~~~~~

If interface doesn't work as expected, or doesn't work at all (!!!) please let us
know. Start by searching our `issue tracker`_ in case someone already reported the
error or behavior you encountered. If you do find an issue, feel free to comment
with details of what you were trying to accomplish and the result. If there isn't
an issue describing your specific situation, open one yourself!

Suggesting New Features
~~~~~~~~~~~~~~~~~~~~~~~

Is there anything you wished interfacea would do? Great, let us know by opening an
issue in our `issue tracker`_! It will be helpful if you write a small description
of the question you want to address, or the problem you are trying to solve. If you
any references or ideas on how to get your idea implemented, let us know as well. We
will be happy to include you as a contributor even if you don't feel comfortable with
writing the code yourself.

Contributing Code
~~~~~~~~~~~~~~~~~

If you want to get your hands dirty and hack away at the core of interfacea, welcome
to the team! To help ourselves stay organized and write cleaner code, we use several
tools such as `git <https://git-scm.com/>`_, `pre-commit <https://pre-commit.com/>`_
and `Github Actions <https://github.com/features/actions>`_.

We make use of branches in git to develop new features and fix bugs in the code. The
:code:`main` branch holds the production version of the code and **should not** ever
be committed to directly. Instead, we create new branches and pull requests to get
new code onto the main branch. If you never used git branches, fear not, there are
plenty of tutorials online, e.g https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow

In short, to start you need either to clone our repository or pull the latest changes
of the :code:`main` branch. Then, you create a new branch and hack away. When you are
done, you should push to your own fork and open a pull request. Feel free to get in
touch with us via an issue or directly via email if you are having trouble!

About writing code itself, as we mentioned before, we make use of a few tools to ease
both your job and ours, when reviewing your contributions. To make things very simple
we use :code:`pre-commit` to handle a series of style checks. To configure it locally
you have to install our developer dependencies. In short, you need to do:

.. code-block:: bash

    # Clone our repository (or your fork)
    git clone https://github.com/joaorodrigues/interfacea.git

    # Create an environment with the main dependencies
    conda env create -f environment.yml
    conda activate interfacea

    # Install the developer dependencies
    conda env update -f devtools/dev-environment.yml

    # Install the pre-commit hooks
    pre-commit install
    pre-commit install-hooks


Style Guides
------------

Coding Style Guide
~~~~~~~~~~~~~~~~~~

The code in interfacea *must* be written in Python 3.x and follow the PEP8. You don't
have to worry about this, :code:`pre-commit` will yell at you if you don't :)

Git(hub) Etiquette
~~~~~~~~~~~~~~~~~~

Please try to follow these guidelines when writing commit messages:
* Write a comprehensive message: no `Fixed bug` messages please!
* Limit the first line to 72 characters.
* Reference issues and pull requests after the first line.


.. _issue tracker: https:/github.com/joaorodrigues/interfacea/issues
