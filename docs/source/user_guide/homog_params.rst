Homogeneous Reaction Parameters
===============================

The parameters of homogeneous/bulk reactions can be provided through
:attr:`echemfem.EchemSolver.homog_params`, a list containing one dictionary for
each reaction. Below is a list of different keys that should appear in each
dictionary

* :Key: ``"stoichiometry"``
  :Type: :py:class:`dict` of :py:class:`int`
  :Description: This entry defines the stoichiometry of the reaction. Each key in this dictionary should be the name of a species (as defined in the ``"name"`` key values of ``conc_params``). The corresponding value is the stoichemetry coefficient for the species in this reaction. It should be an integer: negative for a reactant, positive for a product.
* :Key: ``"forward rate constant"``
  :Type: :py:class:`float` or Firedrake expression
  :Description: Rate constant for the forward reaction.
* :Key: ``"backward rate constant"`` or ``"equilibrium constant"`` (optional)
  :Type: :py:class:`float` or Firedrake expression
  :Description: Rate constant for the backward reaction or equilibrium constant of the reaction. For a reversible reaction, only of one these is required. For an irreversible reaction, leave unset.
* :Key: ``"reference concentration"`` (optional)
  :Type: :py:class:`float` or Firedrake expression
  :Description: Value of 1M in the appropriate units. By setting this, the rate constants are assumed to have the same units (typically 1/s), and the equilibrium constant to be nondimensional.
* :Key: ``"reaction order"`` (optional)
  :Type: :py:class:`dict` of :py:class:`int`
  :Description: This entry can be used to enforce reaction order. Each key in this dictionary should be the name of a species (as defined in the ``"name"`` key values of ``conc_params``). The corresponding value (positive integer) is reaction order of the species in this reaction. If unset, the absolute value of the stoichiometry coeffcient is the reaction order.

Assuming an arbitrary homogeneous reaction system where the forward and
backward rate constants for reaction :math:`i` are :math:`k_i^f` and
:math:`k_i^b`, respectively, and the stoichiometry coefficient of species
:math:`j` in reaction :math:`i` is :math:`s_{j,i}`, then, according to the law
of mass action, the reaction for species :math:`j` is given by

.. math::

   R_j = \sum_i s_{j, i} \left( k_i^f \prod_{n, s_{n, i} < 0} c_n^{-s_{n, i}} - k_i^b \prod_{n, s_{n, i} > 0} c_n^{s_{n, i}} \right),

where :math:`c_n` is the concentration of species :math:`n`. If the equilibrium
constant :math:`K_i` is provided instead of the backward rate constant, then it
is automatically recovered via :math:`k_i^b = \dfrac{k_i^f}{K_i}`. 
The above formula assumes single-step reactions so that the reaction orders are
given by the stoichiometry coefficients. If that is not the case, the exponents
can be replaced by the ``"reaction order"`` inputs.

Note that here the rate constants can have different units, depending on the
number of reactants or products. It is also common to write the reactions in
terms of dimensionless "activities" instead of concentrations. Then, the rate
constants have the same units (typically 1/s) and the equilibrium constants are
dimensionless. With this definition of rate constants, we instead write the
reaction for species :math:`j` as


.. math::

   R_j = \sum_i s_{j, i} c_\mathrm{ref} \left( k_i^f \prod_{n, s_{n, i} < 0} a_n^{-s_{n, i}} - k_i^b \prod_{n, s_{n, i} > 0} a_n^{s_{n, i}} \right),

where :math:`a_n = c_n / c_\mathrm{ref}` is the activity of species :math:`n`
and :math:`c_\mathrm{ref} = 1\text{M}` is a reference concentration, which
needs to be provided in the correct units to ``homog_params``.

Here is an example of a reaction system that can be implement via
``homog_params``. Simplified CO2 - bicarbonate homogeneous reactions:

.. math::

   \mathrm{CO}_2 + \mathrm{OH}^- \xrightleftharpoons[k_1^b]{k_1^f} \mathrm{HCO}_3^-

   \mathrm{HCO}_3^- + \mathrm{OH}^- \xrightleftharpoons[k_2^b]{k_2^f} \mathrm{CO}_3^{2-} + \mathrm{H}_2\mathrm{O}

Here, we assume the dilute solution case where :math:`\mathrm{H}_2\mathrm{O}`
concentration is not tracked. In ``conc_params``, we will have entries for the
other species in this reactions, with ``"name"`` key values: ``"CO2"``,
``"OH"``, ``"HCO3"``, and ``"CO3"``.  We get the following ``homog_params``
list:

.. code::

    [{"stoichiometry": {"CO": -1,
                        "OH": -1,
                        "HCO3": 1},
      "forward rate constant": k1f,
      "backward rate constant": k1b
     },
     {"stoichiometry": {"HCO3": -1,
                        "OH": -1,
                        "CO3": 1},
      "forward rate constant": k2f,
      "backward rate constant": k2b
     }]

Now here is another bicarbonate approximation with a high pH approximation,
where the reaction system is simplified into a single irreversible reaction:

.. math::

   \mathrm{CO}_2 + 2\mathrm{OH}^- \xrightarrow[]{k_1^f} \mathrm{CO}_3^{2-} + \mathrm{H}_2\mathrm{O}

Since this is not actually a single-step reaction, we need to specify the
reaction order for :math:`\mathrm{OH}^-`. We can use the following
``homog_params`` list:

.. code::

    [{"stoichiometry": {"CO2": -1,
                       "OH": -2,
                       "CO3": 1},
     "reaction order": {"OH": 1},
     "forward rate constant": k1f
     }]


As an alternative to the ``homog_params`` interface, the reactions can be
written directly and passed as a function in ``physical_params["bulk
reaction"]`` (:doc:`physical_params`). See ``examples/cylindrical_pore.py`` for
the direct implementation of the high pH approximation of the bicarbonate bulk
reactions. See ``examples/carbonate.py`` and
``examples/carbonate_homog_params.py`` for two equivalent implementations of
the full bicarbonate bulk reactions.
