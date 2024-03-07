Physical Parameters
===================

The general physical parameters need to be provided in
:attr:`echemfem.EchemSolver.physical_params`, a dictionary.

The only key required in all cases is the following one:

:Key: ``"flow"``
:Type: :py:class:`list` of :py:class:`str`
:Description: Each :py:class:`str` represents a transport mechanism. Below is a
              list of tested options and their parameter requirements, and, if
              relevant, the added contribution to the flux of species :math:`k`,
              :math:`\mathbf N_k`.

* ``"diffusion"``: Fickian diffusion

    * ``conc_params``: ``"diffusion coefficient"`` defines :math:`D_k`

    .. math::

        \mathbf N_k \mathrel{+}= -D_k \nabla c_k

* ``"advection"``: Advection

    * ``conc_params``: ``"diffusion coefficient"``
    * Velocity :math:`\mathbf u` defined by overwriting the method :meth:`echemfem.EchemSolver.set_velocity`

    .. math::

        \mathbf N_k \mathrel{+}= c_k \mathbf u

* ``"migration"``: Electromigration (where the mobility constant is given by the Nernst-Einstein relation)

    * ``conc_params``: ``"diffusion coefficient"``, ``"z"``
    * ``physical_params``: ``"F"``, ``"R"``, ``"T"``

    .. math::

        \mathbf N_k \mathrel{+}= -z_k F\frac{D_k}{RT} c_k \nabla \Phi_2.

* ``"electroneutrality"``: Electroneutrality approximation. Implemented through
  a charge-conservation equation as described in `Roy et al, 2022
  <https://doi.org/10.1016/j.jcp.2022.111859>`_. This eliminates a
  concentration from the solution variables. By default, the last concentration
  in ``conc_params`` is eliminated. It can also be specified by passing
  ``"eliminated": True`` in a species' ``conc_params`` entry.

    .. math::

        \sum_k z_k c_k = 0.

* ``"electroneutrality full"``: Electroneutrality approximation, but implemented explicitly.

* ``"porous"``: For flow in a porous media. Diffusion coefficients are replaced
  with the effective diffusion coefficients using the Bruggeman correlation. A
  Poisson equation is solved for the electronic potential :math:`\Phi_1`. If
  using :doc:`echem_params`, the reaction is volumetric and the provided
  reaction is multiplied by the specific surface area, :math:`a_v`.

    * ``physical_params``:

        * ``"solid conductivity"``: defines the electronic conductivity in the bulk of the solid, denoted :math:`\sigma`
        * ``"saturation"``: defines :math:`S`, the liquid saturation. By default :math:`S=1`.
        * ``"porosity"``: defines :math:`\epsilon`
        * ``"specific surface area"``: defines :math:`a_v`.

  .. math::

        D_k^\mathrm{eff} = (\epsilon S)^{1.5} D_k,

        -\nabla \cdot \left((1-\epsilon)^{1.5} \sigma \nabla \Phi_1 \right) = -a_v \sum_j i_j.

* ``"poisson"``: Poisson equation for the ionic potential.

    * ``physical_params``:

        * ``"vacuum permittivity"``: :math:`\epsilon_0`
        * ``"relative permittivity"``: :math:`\epsilon_r`

  .. math::

        -\nabla \cdot \left( \epsilon_0 \epsilon_r \nabla \Phi_2 \right) = F\sum_k z_k c_k,

  neglacting reaction terms. If the vacuum and relative permittivity are not
  defined, the ionic conductivity is used instead :math:`\kappa = \sum_k z_k^2
  F^2 \frac{D_k}{RT} c_k`.

* ``"finite size"``: Finite-size effect for ion interaction. Choosing
  ``"poisson"``, ``"diffusion"``, ``"migration"``, and ``"finite size"`` gives
  the Generalized Poisson-Nernst-Planck model (GMPNP).

    * ``conc_params``: ``"solvated diameter"``, :math:`a_k` of species :math:`k`
    * ``physical_params``: ``"Avogadro constant"``, :math:`N_A`

    .. math::

        \mathbf N_k \mathrel{+}= -D_k c_k \left(\frac{N_A \sum_j a_j^3 \nabla
        c_j}{1- N_A \sum_j a_j^3 c_j}\right)

Below is a list of other keys that can appear in the dictionary


* :Key: ``"bulk reaction"``
  :Type: a function
  :Description: Homogeneous/bulk reactions to be added to the right-hand side of mass conservation equations. These can instead be set using ``homog_params`` (:doc:`homog_params`).

    Args:
        u: solution state. The value of the different concentrations can be recovered through ``u([self.i_c["species name"])`` within a ``echemfem.EchemSolver`` object.

    Returns:
        List of length equal to the number of species, each entry being the reaction term of the corresponding species. The order of the reaction terms must be the same as the order of the species in ``conc_params``, which can also be found through ``self.i_c``.

* :Key: ``"F"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Faraday Constant
* :Key: ``"R"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Molar Gas Constant
* :Key: ``"T"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Absolute Temperature
* :Key: ``"porosity"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Void fraction of the porous medium. Used when ``"porous"`` is in ``physical_params["flow"]``.
* :Key: ``"solid conductivity"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Electronic conductivity of the bulk solid material. Used when ``"porous"`` is in ``physical_params["flow"]``.
* :Key: ``"specific surface area"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Specific surface area of the porous medium. Used when ``"porous"`` is in ``physical_params["flow"]``, if :doc:`echem_params` are used.
* :Key: ``"saturation"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Ratio of the void fraction occupied by liquid. Default value of ``1.0``. Used when ``"porous"`` is in ``physical_params["flow"]``.
* :Key: ``"vacuum permittivity"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Used when ``"poisson"`` is in ``physical_params["flow"]``.
* :Key: ``"relative permittivity"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Used when ``"poisson"`` is in ``physical_params["flow"]``.
* :Key: ``"Avogadro constant"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Used when ``"finite size"`` is in ``physical_params["flow"]``.
* :Key: ``"U_app"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Applied potential. Used for ``"applied"``, ``"liquid applied"``, and ``"robin"`` :doc:`boundary_conditions`.
* :Key: ``"gap capacitance"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Used for ``"robin"`` :doc:`boundary_conditions`.
* :Key: ``"surface charge density"``
  :Type: :py:class:`float`, firedrake expression
  :Description: Used for ``"poisson neumann"`` :doc:`boundary_conditions`.

