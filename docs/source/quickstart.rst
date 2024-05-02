Quickstart Guide
================

.. _installation:

Installation
------------

First, please install the open-source finite element library `Firedrake <https://www.firedrakeproject.org>`_.
On a Mac with Homebrew installed or on an Ubuntu workstation with sudo access, it can be installed using the default configuration as detailed `here <https://www.firedrakeproject.org/download.html>`_.
To get started quickly, a `Docker image <https://hub.docker.com/r/firedrakeproject/firedrake>`_ is also available.
Note that you may need to increase the default memory allocation for Docker in order to run some of the EchemFEM examples.

EchemFEM is hosted on `GitHub <https://github.com/LLNL/echemfem>`_, and should be cloned from there.

To use EchemFEM, first install it using pip within the Firedrake virtual environment. From the echemfem parent directory:

.. code-block:: console

   (firedrake) $ pip install -e .

To test your installation, you can run the tests in `echemfem/tests <https://github.com/LLNL/echemfem/tree/main/tests>`_ using pytest.
This can take some time, so it is probably enough to run this test:

.. code-block:: console

   (firedrake) $ pytest test_advection_diffusion_migration.py

Alternatively, you can run some of the examples below.

Running Examples
----------------

To get started, a demo for a flow reactor is available :doc:`here <demos/index>`.
Several examples can be found in `echemfem/examples <https://github.com/LLNL/echemfem/tree/main/examples>`_.
Examples using the generalized modified Poisson-Nernst-Planck model (GMPNP) can be found in `FireCat <https://github.com/LLNL/firecat>`_, where it is coupled with a microkinetics model.

To run an example:

.. code-block:: console

   (firedrake) $ python example_name.py

And to run it in parallel with ``n_proc`` processors:

.. code-block:: console

   (firedrake) $ mpiexec -n n_proc python example_name.py

Visualization
-------------

The solution fields are stored for visualization with `Paraview <https://www.paraview.org>`_ in ``results/collection.pvd``.
Alternatively, 1D and 2D functions can be visualized using `matplotlib <https://matplotlib.org>`_. 
More about `Firedrake visualization <https://www.firedrakeproject.org/visualisation.html>`_.
