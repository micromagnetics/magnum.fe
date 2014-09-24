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

For a detailed discussion on the LLG see [1]. magnum.fe implements different algorithms for the integration of the LLG. Note that some algorithms depend on a specific linear backend and/or additional packages.

References:
  [1] Brown Jr, W. F. Micromagnetics, 1963. Interscience, New York.

:class:`LLGAlougesProject`
++++++++++++++++++++++++++

.. autoclass:: LLGAlougesProject
   :members:

:class:`LLGAlougesLagrange`
+++++++++++++++++++++++++++

.. autoclass:: LLGAlougesLagrange
   :members:
