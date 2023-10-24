Boundary Conditions
===================

The boundary conditions need to set through
:meth:`echemfem.EchemSolver.set_boundary_markers`, which sets a dictionary
containing boundary condition names and their corresponding boundary id.

In the equations below, we have:

* :math:`c_k`, the concentration of species :math:`k`
* :math:`\Phi_2`, ionic potential
* :math:`\Phi_1`, electronic potential
* :math:`\mathbf N_k`, the flux of species :math:`k`
* :math:`c_{k,\mathrm{bulk}}`, the bulk concentration of species :math:`k`
* :math:`\mathbf n`, the unit outward normal vector
* :math:`\mathbf u`, the velocity

Here are the different options for ``boundary_markers`` keys:

* ``"inlet"``: Inlet boundary condition defined as

.. math::

    \mathbf N_k \cdot \mathbf n = c_{k,\mathrm{bulk}} \mathbf u \cdot \mathbf n.

* ``"outlet"``: Outlet boundary condition defined as

.. math::

    \mathbf N_k \cdot \mathbf n = c_k \mathbf u \cdot \mathbf n.

* ``"bulk dirichlet"``: Dirichlet boundary condition for concentrations using the value ``"bulk"`` provided in :doc:`conc_params` such that

.. math::

    c_k = c_{k,\mathrm{bulk}}.

* ``"bulk"``: Robin boundary condition for concentrations, where :math:`K_{k,MT}` is the value ``"mass transfer coefficient"`` if provided in :doc:`conc_params` for species :math:`k`. Also homogeneous Dirichlet boundary condition of the ionic potential, as follows

.. math::

    \begin{align}
    \mathbf N_k \cdot \mathbf n &= K_{k,MT} (c_k - c_{k,\mathrm{bulk}}), \\
    \Phi_2 &= 0.
    \end{align}

* ``"neumann"``: Neumann boundary condition defined through the method :meth:`echemfem.EchemSolver.neumann`.
* ``"gas"``:  Dirichlet boundary condition for concentrations using the value ``"gas"`` if provided in :doc:`conc_params` for species :math:`k`.
* ``"applied"``: Dirichlet boundary condition for ionic potential, if this is a non-porous case, and for the electronic potential, if this is a porous case. It uses the value ``"U_app"`` given in :doc:`physical_params`.
* ``"liquid applied"``: Dirichlet boundary condition for ionic potential using the value ``"U_app"`` given in :doc:`physical_params`. This is the same as ``"applied"`` in the non-porous case.
* ``"robin"``: Robin boundary condition for the ionic potential when using the Poisson equation, typically used for GMPNP. We pass :math:`U - U^\mathrm{PZC}` together as ``"U_app"`` in :doc:`physical_params`, as well as ``"gap capacitance"`` for :math:`C_\mathrm{gap}`.

.. math::
   \epsilon_0\epsilon_\mathrm{r} \nabla \Phi_2 \cdot \mathbf{n}= C_\mathrm{gap}\left(U - U^\mathrm{PZC} - \Phi_2\right).

* ``"poisson neumann"``: Neumann boundary condition for the ionic potential when using the Poisson equation, where :math:`\sigma` is the value ``"surface charge density"`` in :doc:`physical_params`.

.. math::
   \epsilon_0\epsilon_\mathrm{r} \nabla \Phi_2 \cdot \mathbf{n}= \sigma

* Charge-transfer reactions: a custom :py:class:`str` can be passed to name the electrodes used in surface charge-transfer reactions defined using :doc:`echem_params`.

