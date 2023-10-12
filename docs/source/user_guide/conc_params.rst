Concentration Parameters
========================

The physical parameters of each species need to be provided
:attr:`echemfem.EchemSolver.conc_params`, a list containing one dictionary for
each species. Below is a list of different keys that can appear in each dictionary

* :Key: ``"name"``
  :Type: ``str``
  :Description: Species name. E.g. ``"CO2"``
* :Key: ``"bulk"``
  :Type: ``float``, firedrake expression
  :Description: Concentration at the "bulk". This value is used in :meth:`echemfem.EchemSolver.setup_solver`


