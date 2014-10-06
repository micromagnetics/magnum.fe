.. module:: magnumfe
.. _spin-diffusion:

Spin Diffusion
==============

The interaction of spin polarized currents with the magnetization can be described by a spin-accumulation model described in [Zhang2002]_ and [Garcia2007]_.
In this model the spin accumulation :math:`\vec{s}(\vec{x})` exerts a torque on the magnetization by an additional term to the LLG

.. math::
  \frac{\partial \vec{m}}{\partial t} =
  - \gamma \vec{m} \times (\vec{h}_\text{eff} + \frac{c}{\mu_0} \vec{s})
  + \alpha \vec{m} \times \frac{\partial \vec{m}}{\partial t}

where :math:`c` is a coupling constant. The spin accumulation itself depends on the magnetization by the diffusion equation

.. math::
  \frac{\partial \vec{s}}{\partial t} =
  - \vec{\nabla} \cdot \mat{J}_\text{s}
  - 2 D_0 \left[
    \frac{\vec{s}}{\lambda_\text{sf}^2}
    + \frac{\vec{s} \times \vec{m}}{\lambda_\text{J}^2}
  \right]

where :math:`D_0` is the diffusion constant, :math:`\lambda_\text{sf}` is the characteristic length for spin-flip relaxation, :math:`\lambda_\text{J}` depends on the electron's mean free path and :math:`\mat{J}_\text{s}` is the matrix valued spin current given by

.. math::
  \mat{J}_\text{s} = 
  \frac{\beta \mu_\text{B}}{e} \vec{m} \otimes \vec{J}_\text{e}
  - 2 D_0 \left[
    \vec{\nabla} \vec{s}
    - \beta \beta' \vec{m} \otimes \left( (\vec{\nabla}\vec{s})^T \vec{m} \right)
  \right].

Here :math:`\beta` and :math:`beta'` are dimensionless polarization parameters and :math:`\mat{J}_\text{e}` is the electrical current density. As shown in [Abert2014b]_ this method is able to reproduce the results of both the spin torque model by Slonczewski, see [Slonczewski1996]_, and the model by Zhang and Li, see [Zhang2004]_.

:class:`SpinDiffusion`
++++++++++++++++++++++

.. autoclass:: SpinDiffusion
   :members:


.. [Zhang2002] Zhang, S., Levy, P. M., & Fert, A. (2002). Mechanisms of spin-polarized current-driven magnetization switching. Physical review letters, 88(23), 236601.
.. [Garcia2007] García-Cervera, C. J., & Wang, X. P. (2007). Spin-polarized currents in ferromagnetic multilayers. Journal of Computational Physics, 224(2), 699-711.
.. [Abert2014b] Abert, C., Ruggeri M., Bruckner F., Vogler C., Hrkac G., Praetorius D., & Suess D. (2014). Self-consistent micromagnetic simulations including spin diffusion effects. In preparation
.. [Slonczewski1996] Slonczewski, J. C. (1996). Current-driven excitation of magnetic multilayers. Journal of Magnetism and Magnetic Materials, 159(1), L1-L7.
.. [Zhang2004] Zhang, S., & Li, Z. (2004). Roles of nonequilibrium conduction electrons on the magnetization dynamics of ferromagnets. Physical Review Letters, 93(12), 127204.
.. [Abert2014a] C. Abert, G. Hrkac, M. Page, D. Praetorius, M. Ruggeri, & D. Suess.  Spin-polarized transport in ferromagnetic multilayers: An unconditionally convergent FEM integrator.  Comput. Math. Appl., 68, 6, 639-654, 2014. 
