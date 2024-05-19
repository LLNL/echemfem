Electrochemical Parameters
==========================

The parameters of each charge-transfer reaction can be provided through
:attr:`echemfem.EchemSolver.echem_params`, a list containing one dictionary for
each charge-transfer reaction. Below is a list of different keys that can
appear in each dictionary

* :Key: ``"reaction"``
  :Type: a function
  :Description: the current density from a charge-transfer reaction

    Args:
        u: solution state. The value of the different concentrations can be recovered through ``u([self.i_c["species name"])`` within a ``echemfem.EchemSolver`` object.

    Returns:
        Current density :math:`i_k`.

* :Key: ``"electrons"``
  :Type: :py:class:`float`, :py:class:`int`, firedrake expression
  :Description: the number of electrons :math:`n_k` being transferred in the reaction.
* :Key: ``"stoichiometry"``
  :Type: :py:class:`dict`
  :Description: This entry defines the stoichiometry of the reaction. Each key in this dictionary should be the name of a species (as defined in the ``"name"`` key values of ``conc_params``). The corresponding value is the stoichiometry coefficient :math:`s_{j,k}` for the species in this reaction. It should be an integer: negative for a reactant, positive for a product.
* :Key: ``"boundary"``
  :Type: :py:class:`str`
  :Description: In the non-porous case, this entry provides the name of the boundary where the charge-transfer reaction happens. This should correspond to a key in the dictionary ``self.boundary_markers`` containing the :doc:`boundary_conditions`.

Here is an example of a reaction that can be implement via ``echem_params``. CO reduction for ethanol:

.. math::

   2\mathrm{CO} + 7 \mathrm{H}_2 \mathrm{O} + 8 \mathrm{e}^- \rightarrow \mathrm{CH}_3\mathrm{CH}_2\mathrm{OH}_{(\mathrm{aq})} + 8 \mathrm{OH}^-

Here, we assume the dilute solution case where :math:`\mathrm{H}_2\mathrm{O}`
concentration is not tracked. In ``conc_params``, we will have entries for the other species in this reactions, with ``"name"`` key values: ``"CO"``, ``"OH"``, and ``"C2H6O"``. Assuming we also defined a function ``reaction_C2H6O`` for the current density, we get the following ``echem_params`` entry:

.. code::

    {"reaction": reaction_C2H6O,
     "electrons": 8,
     "stoichiometry": {"CO": -2,
                       "C2H6O": 1,
                       "OH": 8},
     }

In the non-porous case, the reactions are surface reactions, which result in the following boundary conditions for the ionic current :math:`\mathbf{i}_2`, and mass flux :math:`\mathbf{N}_j` of species :math:`j`

.. math::

    \mathbf{i}_2 \cdot \mathbf{n} = -\sum_k i_k,

    \mathbf{N}_j \cdot \mathbf n = -\sum_k \frac{s_{j,k} i_k}{n_k F},

where :math:`k` represent each charge-transfer reaction in the system, and :math:`F` is the Faraday constant, which must be passed through :doc:`physical_params`.

In the porous case, i.e. when ``"porous"`` is passed through ``"flow"`` in :doc:`physical_params`, the reactions are volumetric reactions, which are added to the right-hand sides of the equations. In the absence of homogeneous bulk reactions, we have

.. math::

    \nabla \cdot \mathbf{i}_1 = -\nabla \cdot \mathbf{i}_2 = - a_v \sum_k i_k,

    \nabla \cdot \mathbf{N}_j =  a_v \sum_k \frac{s_{j,k} i_k}{n_k F},

where :math:`a_v` is passed as ``"specific surface area"`` through :doc:`physical_params`.
