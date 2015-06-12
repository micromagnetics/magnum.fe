.. module:: magnumfe

Landau-Lifshitz-Gilbert Equation
================================

The Landau-Lifshitz-Gilbert equation (LLG) describes the motion of the magnetization configuration in micromagnetics:

.. math::
  \partial_t \vec{m} =
  - \gamma ( \vec{m} \times \vec{H}_\text{eff} )
  + \alpha ( \vec{m} \times \partial_t \vec{m} )

where the effective field :math:`\vec{H}_\text{eff}` is given by the variational derivative of the total free energy

.. math::
  \vec{H}_\text{eff} =
  - \frac{1}{\mu_0 M_\text{s}} \frac{\delta U(\vec{m})}{\delta \vec{m}}.

For a detailed discussion on the LLG see [Brown1963]_. magnum.fe implements different algorithms for the integration of the LLG. Note that some algorithms depend on a specific linear backend and/or additional packages.

:class:`LLGAlougesProject`
++++++++++++++++++++++++++

.. autoclass:: LLGAlougesProject
   :members:

:class:`LLGAlougesLagrange`
+++++++++++++++++++++++++++

.. autoclass:: LLGAlougesLagrange
   :members:


.. [Brown1963] Brown Jr, W. F. Micromagnetics, 1963. Interscience, New York.
.. [Suess2002] Suess, D., Tsiantos, V., Schrefl, T., Fidler, J., Scholz, W., Forster, H., ... & Miles, J. J. (2002). Time resolved micromagnetics using a preconditioned time integration method. Journal of Magnetism and Magnetic Materials, 248(2), 298-311.
.. [Alouges2008] Alouges, F. (2008). A new finite element scheme for Landau-Lifchitz equations. Discrete Contin. Dyn. Syst. Ser. S, 1(2), 187-196.
.. [Abert2013b] Abert, C., Exl, L., Bruckner, F., Drews, A., & Suess, D. (2013). magnum. fe: A micromagnetic finite-element simulation code based on FEniCS. Journal of Magnetism and Magnetic Materials, 345, 29-35.
.. [Goldenits2012] Goldenits, P., Hrkac, G., Praetorius, D., & Suess, D. (2012, February). An effective integrator for the Landau-Lifshitz-Gilbert equation. In Proceedings of Mathmod 2012 Conference.
