User Guide
==========

This guide will help a user understand existing examples, and design new ones.

The first step of writing an EchemFEM script is to import the ``EchemSolver`` class:

.. code-block::

   from echemfem import EchemSolver

or for a transient simulation, the ``TransientEchemSolver`` class:

.. code-block::

   from echemfem import TransientEchemSolver

These are abstract classes, which we will use as the base class for a specific model.
To create a customized solver, the user should set the following inputs:

.. toctree::
   :maxdepth: 2

   conc_params
   physical_params
   boundary_conditions
   echem_params
   homog_params

Here is a barebone example of how a user might define their own solver.
Several concrete examples can be found in `echemfem/examples
<https://github.com/LLNL/echemfem/tree/main/examples>`_.

.. code-block::

    class MySolver(EchemSolver): # or TransientEchemSolver

        def __init__(self):
            # Here define all custom parameters that may require attributes

            super().__init__(...) # with appropriate arguments

        def set_boundary_markers(self):
            self.boundary_markers = ...

        # and some other methods that need to be defined

Then, to run the simulation, create the object and run the ``solve`` method.

.. code-block::

   solver = MySolver()
   solver.solve()

For transient cases, a temporal grid defined as a  ``list`` or
``numpy.ndarray`` must be provided. For example,

.. code-block::

   import numpy as np
   times = np.linspace(0, 11, 1)  # 10 timesteps of size 0.1
   solver.solve(times)

To generate non-trivial velocity fields, some fluid flow solvers are available:

.. toctree::
   :maxdepth: 2

   fluid_solver

