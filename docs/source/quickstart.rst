Quickstart Guide
================

.. _installation:

Installation
------------

First, please install the open-source finite element library `Firedrake <https://www.firedrakeproject.org/download.html>`_.

EchemFEM is hosted on `GitHub <https://github.com/LLNL/echemfem>`_, and should be cloned from there.

To use EchemFEM, first install it using pip within the Firedrake virtual environment. From the echemfem parent directory:

.. code-block:: console

   (firedrake) $ pip install -e .

To test your installation, you can run the tests in `echemfem/tests <https://github.com/LLNL/echemfem/examples>`_ using pytest.
This can take some time, so it is probably enough to run this test:

.. code-block:: console

   (firedrake) $ pytest test_advection_diffusion_migration.py

Alternatively, you can run some of the examples below.

Running Examples
----------------

Several examples can be found in `echemfem/examples <https://github.com/LLNL/echemfem/examples>`_.
Examples using the generalized modified Poisson-Nernst-Planck model (GMPNP) can be found in [FireCat](https://github.com/LLNL/firecat), where it is coupled with a microkinetics model.
