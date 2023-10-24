Concentration Parameters
========================

The physical parameters of each species need to be provided in
:attr:`echemfem.EchemSolver.conc_params`, a list containing one dictionary for
each species. Below is a list of different keys that can appear in each dictionary.
Only the first key ``"name"`` is required for every case.

* :Key: ``"name"``
  :Type: :py:class:`str`
  :Description: Species name. E.g. ``"CO2"``
  :Uses: * Name of field in pvd output
         * Used to get index in the solution vector.
* :Key: ``"bulk"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Concentration at the "bulk".
  :Uses: * Initial guess for the concentrations in :meth:`echemfem.EchemSolver.setup_solver`
         * Concentration value for ``"bulk dirichlet"``, ``"bulk"``, and ``"inlet"`` :doc:`boundary_conditions`
* :Key: ``"gas"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Concentration at the gas interface for the ``"gas"`` :doc:`boundary_conditions`.
* :Key: ``"diffusion coefficient"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Diffusion coeffcient of the species.
* :Key: ``"z"``
  :Type: :py:class:`float`, :py:class:`int`, firedrake expression
  :Description: Charge number of the species.
* :Key: ``"mass transfer coefficient"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Coefficient used in the ``"bulk"`` :doc:`boundary_conditions`.
* :Key: ``"solvated diameter"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Solvated diameter of the ionic species if using ``"finite size"`` in :doc:`physical_params`.
* :Key: ``"eliminated"``
  :Type: :py:class:`bool`
  :Description: The species to be eliminated via the electroneutrality approximation if using ``"electroneutrality"`` in :doc:`physical_params`.
* :Key: ``"C_ND"``
  :Type: :py:class:`float`
  :Description: For the nondimensionalization used in `Roy et al, 2022 <https://doi.org/10.1016/j.jcp.2022.111859>`_. These weights are used to get the nondimensional charge conservation equation. In the paper :math:`C_{ND} = c_k^\mathrm{in} / c_\mathrm{ref}`.

