Fluid Flow Solvers
===================

For convenience, we provide a simple implementation of fluid flow equations in
:class:`echemfem.FlowSolver`, which can be used to obtain a velocity field for
the :class:`echemfem.EchemSolver`.

Navier-Stokes Solver
--------------------

The :class:`echemfem.NavierStokesFlowSolver` class contains an implementation
for the incompressible Navier-Stokes equations. Both nondimensional and
dimensional equations are available. Physical parameters are passed as a
:py:class:`dict` in the ``fluid_params`` argument.

To use the nondimensional version, the user must pass the following key:

* ``"Reynolds number"``

For the dimensional version, the user must pass the following keys:

* ``"density"``

* ``"dynamic viscosity"`` or ``"kinematic viscosity"``

Navier-Stokes-Brinkman Solver
-----------------------------

The :class:`echemfem.NavierStokesBrinkmanFlowSolver` class contains an
implementation for the incompressible Navier-Stokes-Brinkman equations. Both
nondimensional and dimensional equations are available. In addition to the
parameters provided for Navier-Stokes, the physical parameters below must be
provided.

To use the nondimensional version, the user must also pass the following key:

* ``"Darcy number"``

For the dimensional version, the user must also pass the following keys:

* ``"permeability"`` or ``"inverse permeability"``

* Optional: ``"effective kinematic viscosity"``

Boundary conditions
-------------------

A :py:class:`dict` containing boundary markers is passed when creating a
``FlowSolver`` object. Below are the different options for ``boundary_markers``
keys. The velocity is denoted as :math:`\mathbf u` and pressure as :math:`p`.

* ``"no slip"``: Sets velocity to zero, :math:`\mathbf u = 0`.

* ``"inlet velocity"``: Sets inlet velocity, :math:`\mathbf u = \mathbf
  u_\mathrm{in}`, which is passed in ``fluid_params`` as ``"inlet velocity"``.

* ``"outlet velocity"``: Sets outlet velocity, :math:`\mathbf u = \mathbf
  u_\mathrm{out}`, which is passed in ``fluid_params`` as ``"outlet
  velocity"``.

* ``"inlet pressure"``: Sets inlet velocity, :math:`p = p_\mathrm{in}`, which
  is passed in ``fluid_params`` as ``"inlet pressure"``.

* ``"outlet pressure"``: Sets outlet pressure, :math:`p = p_\mathrm{out}`,
  which is passed in ``fluid_params`` as ``"outlet pressure"``.
