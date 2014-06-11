.. module:: magnumfe

Landau-Lifshitz-Gilbert Equation
================================

The Landau-Lifshitz-Gilbert equation (LLG) describes the motion of the magnetization configuration in micromagnetics:

.. math::
  \partial_t \boldsymbol{m} =
  - \gamma ( \boldsymbol{m} \times \boldsymbol{H}_\text{eff} )
  + \alpha ( \boldsymbol{m} \times \partial_t \boldsymbol{m} )

where the effective field :math:`\boldsymbol{H}_\text{eff}` is given by the variational derivative of the total free energy

.. math::
  \boldsymbol{H}_\text{eff} =
  - \frac{1}{\mu_0 M_\text{s}} \frac{\delta U(\boldsymbol{m})}{\delta \boldsymbol{m}}.

For a detailed discussion on the LLG see [1]. magnum.fe implements different algorithms for the integration of the LLG. Note that some algorithms depend on a specific linear backend and/or additional packages.

[1] Brown Jr, W. F. Micromagnetics, 1963. Interscience, New York.

:class:`LLGAlougesProject`
++++++++++++++++++++++++++

.. autoclass:: LLGAlougesProject
   :members:

:class:`LLGAlougesLagrange`
+++++++++++++++++++++++++++

.. autoclass:: LLGAlougesLagrange
   :members:
