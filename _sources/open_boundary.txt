.. module:: magnumfe
.. _open-boundary:

Open Boundary Problems
======================

Multiple effective-field terms, namely the demagnetization field (:class:`DemagField`) and the Oersted field (:class:`OerstedField`) involve the solution of so called open-boundary problems of the form

.. math::
  \Delta u                                         &= \nabla \cdot \vec{f}     \quad \text{in} \quad \omega \\
  \Delta u                                         &= 0                        \quad \text{in} \quad \mathbb{R}^3 \setminus \omega \\
  [u]                                              &= 0                        \quad \text{on} \quad \partial\omega \\
  \left[\frac{\partial u}{\partial \vec{n}}\right] &= - \vec{f} \cdot \vec{n}  \quad \text{on} \quad \partial\omega \\
  u(\vec{x})                                       &= \mathcal{O}(1/|\vec{x}|) \quad \text{if} \quad |\vec{x}| \rightarrow \infty

different algorithms are available for the solution of these problems.

.. _shell-transformation:

Shell-Transformation Method
+++++++++++++++++++++++++++
  
The shell-transformation method for the solution of open-boundary problems with the finite-element method was introduced in [1].
The region of interest :math:`\omega{notrans}` is surrounded by a finite shell :math:`\omega{trans}`.
A bijective transformation is constructed that maps the shell region :math:`\omega{trans}` onto the complete exterior region :math:`\mathbb{R}^3 \setminus \omega_\text{notrans}`.
The open-boundary conditions can thus be implemented by requiring homogeneous Dirichlet boundary condition of the outer boundary of the shell :math:`\omega{trans}`, resulting in the following weak formulation

.. math::
  \int_{\omega_\text{notrans}} \nabla u \cdot \nabla v \dx + \int_{\omega_\text{trans}} (\nabla u)^T \mat{g} \nabla v \dx = \int_{\omega} \vec{f} \cdot \nabla v \dx

with the metric tensor :math:`\mat{g}` defined as

.. math::
  \mat{g} = (\mat{J}^{-1})^T | \det \mat{J} | \mat{J}^{-1}

where :math:`\mat{J}` denotes the Jacobian of the transformation.
The untransformed region :math:`\omega_\text{notrans}` does not have to coincide with the sample region :math:`\omega` of the original problem, but it holds :math:`\omega \subseteq \omega_\text{notrans}`.

magnum.fe uses cuboid shells, i.e. `\omega_\text{notrans}` has to be of cuboidal shape. The applied method is described in detail in [2].

.. _fem-bem:

Hybrid FEM-BEM Method (Fredkin and Koehler)
+++++++++++++++++++++++++++++++++++++++++++

A hybrid FEM-BEM approach for the solution of open-boundary problems was propsed by Fredkin and Koehler [3]. Consider the following splitting of the solution :math:`u`:

.. math::
    u = u_1 + u_2

:math:`u_1` is defined by

.. math::
  \Delta u_1                            &= \nabla \cdot \vec{f}    \quad \text{in} \quad \omega \\
  \frac{\partial u_1}{\partial \vec{n}} &= - \vec{f} \cdot \vec{n} \quad \text{on} \quad \partial\omega \\
  u_1                                   &= 0                       \quad \text{in} \quad \mathbb{R}^3 \setminus \omega.

This Neuman problem within :math:`\omega` is solved with the finite-element method.
While :math:`u_1` solves for the right-hand side :math:`\nabla \cdot \vec{f}` and fullfills the jump condition of the normal derivative :math:`- \vec{f} \cdot \vec{n}` it is not continuous across :math:`\partial \omega`.
This jump is compensated by :math:`u_2` which is defined as

.. math::
  \Delta u_2                                         &= 0                        \quad \text{in} \quad \omega \\
  [u_2]                                              &= - [u_1]                  \quad \text{on} \quad \partial \omega \\
  \left[\frac{\partial u_2}{\partial \vec{n}}\right] &= 0                        \quad \text{on} \quad \partial\omega \\
  u_2(\vec{x})                                       &= \mathcal{O}(1/|\vec{x}|) \quad \text{if} \quad |\vec{x}| \rightarrow \infty

This system is solved by the double-layer potential

.. math::
  \begin{align}
    u_2 = \int_{\partial \omega} u_1 \frac{\partial}{\partial \vec{n}} \frac{1}{|\vec{x} - \vec{x}'|} \dx
  \end{align}

The double-layer potential is computed on the boundary :math:`\partial \omega` with the boundary-element method.
These values are used as Dirichlet boundary conditions to solve :math:`u_2` within :math:`\omega` with the finite-element method.
For a detailed description of the method, see [3].

References:
  [1] Brunotte, X., Meunier, G., & Imhoff, J. F. (1992). Finite element modeling of unbounded problems using transformations: a rigorous, powerful and easy solution. IEEE Transactions on Magnetics, 28(2), 1663-1666.

  [2] Abert, C., Exl, L., Selke, G., Drews, A., & Schrefl, T. (2013). Numerical methods for the stray-field calculation: A comparison of recently developed algorithms. Journal of Magnetism and Magnetic Materials, 326, 176-185.

  [3] Fredkin, D. R., & Koehler, T. R. (1990). Hybrid method for computing demagnetizing fields. Magnetics, IEEE Transactions on, 26(2), 415-417.
