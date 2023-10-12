User Guide
==========

This guide will help a user understand existing examples, and design new ones.

The first step of writing an EchemFEM script is to import the ``EchemSolver`` class:

.. code-block::

   from echemfem import EchemSolver

This is an abstract class, which we will use as the base class for a specific model.
To create a customized solver, the user should set the following inputs:

.. toctree::
   :maxdepth: 2

   conc_params
   physical_params
   boundary_conditions
   echem_params

