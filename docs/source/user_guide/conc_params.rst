Concentration Parameters
========================

The physical parameters of each species need to be provided
:attr:`echemfem.EchemSolver.conc_params`, a list containing one dictionary for
each species. Below is a list of different keys that can appear in each dictionary

* :Key: ``"name"``
  :Type: :py:class:`str`
  :Description: Species name. E.g. ``"CO2"``
  :Uses: Name of field in pvd output. Used to get index in the solution vector.
* :Key: ``"bulk"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Concentration at the "bulk". This value is used in :meth:`echemfem.EchemSolver.setup_solver`

.. note::

   This page is under construction.
