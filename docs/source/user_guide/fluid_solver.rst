Fluid Flow Solvers
===================

For convenience, we provide a simple implementation of fluid flow equations in
:class:`echemfem.FlowSolver`, which can be used to obtain a velocity field for
the :class:`echemfem.EchemSolver`.

Navier-Stokes Solver
--------------------

The :class:`echemfem.NavierStokesFlowSolver` class contains an implementation
for the steady-state incompressible Navier-Stokes equations. Both
nondimensional and dimensional equations are available.

We are solving for velocity :math:`\mathbf u` and pressure :math:`p`. The
dimensional version of the momentum equation is

.. math::

   -\nu \nabla^2 \mathbf u + \mathbf u \cdot \nabla \mathbf u + \frac{1}{\rho} \nabla p = 0,

where :math:`\nu` is the kinematic viscosity and :math:`\rho` is the density.
It is also common to use the dynamic viscosity :math:`\mu = \nu \rho`

The nondimensional version is

.. math::

   - \frac{1}{\mathrm{Re}}\nabla^2 \mathbf u + \mathbf u \cdot \nabla \mathbf u + \nabla p = 0,

where all quantities have been nondimensionalized, and :math:`\mathrm{Re}` is the Reynolds number. In both cases, we also have the incompressibility condition

.. math::

   \nabla \cdot \mathbf{u} = 0.

Physical parameters are passed as a
:py:class:`dict` in the ``fluid_params`` argument.

For the dimensional version, the user must pass the following keys:

* ``"density"``

* ``"dynamic viscosity"`` or ``"kinematic viscosity"``

To use the nondimensional version, the user must pass the following key:

* ``"Reynolds number"``

Navier-Stokes-Brinkman Solver
-----------------------------

The :class:`echemfem.NavierStokesBrinkmanFlowSolver` class contains an
implementation for the steady-state incompressible Navier-Stokes-Brinkman
equations. Both nondimensional and dimensional equations are available.

The dimensional version of the momentum equation is

.. math::

   -\nabla\cdot\left(\nu_\mathrm{eff} \nabla\mathbf u \right)+ \mathbf u \cdot \nabla \mathbf u + \frac{1}{\rho} \nabla p + \nu K^{-1} \mathbf u = 0,

where :math:`\nu_\mathrm{eff}` is the effective viscosity in the porous medium,
and :math:`K`, its permeability. The inverse permeability can be provided
directly in cases where it is zero in some regions, i.e. liquid-only regions.
It is common to take :math:`\nu_\mathrm{eff}=\nu` for simplicity (the default
here).

The nondimensional implementation currently assume :math:`\nu_\mathrm{eff}=\nu`
and :math:`K>0`, such that


.. math::

   - \frac{1}{\mathrm{Re}}\nabla^2 \mathbf u + \mathbf u \cdot \nabla \mathbf u + \nabla p + \frac{1}{\mathrm{Re}\mathrm{Da}} \mathbf u = 0,

where all quantities have been nondimensionalized and :math:`\mathrm{Da}` is the Darcy number.

In addition to the parameters provided for Navier-Stokes, the physical
parameters below must be provided.

For the dimensional version, the user must also pass the following keys:

* ``"permeability"`` or ``"inverse permeability"``

* Optional: ``"effective kinematic viscosity"`` or ``"effective dynamic viscosity"``

To use the nondimensional version, the user must also pass the following key:

* ``"Darcy number"``

Boundary conditions
-------------------

A :py:class:`dict` containing boundary markers is passed when creating a
``FlowSolver`` object. Below are the different options for ``boundary_markers``
keys.

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
