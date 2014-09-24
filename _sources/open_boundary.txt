.. module:: magnumfe

Open Boundary Problems
======================

Multiple effective-field terms, namely the demagnetization field (:class:`DemagField`) and the Oersted field (:class:`OerstedField`) involve the solution of so called open-boundary problems of the form

.. math::
  \begin{align}
    \Delta u     &= f(\vec{x})             &\text{in}&\quad \omega \\
    \Delta u     &= 0                      &\text{in}&\quad \mathbb{R}^3 \setminus \omega \\
    \partial_n u &= g(\vec{x})             &\forall  &\quad \vec{x} \in \partial \omega \\
    u(\vec{x}) &= \mathcal{O}(1/|\vec{x}|) &\text{if}&\quad |\vec{x}| \rightarrow \infty
  \end{align}

different algorithms are available for the solution of these problems.

Shell Transformation
++++++++++++++++++++


Hybrid FEM-BEM Method (Fredkin and Koehler)
+++++++++++++++++++++++++++++++++++++++++++

A hybrid FEM-BEM approach for the solution of open-boundary problems was propsed by Fredkin and Koehler [1]. Consider the following splitting of the solution :math:`u`:

.. math::
    u = u_1 + u_2

:math:`u_1` is zero within :math:`\omega` and solves the Neumann problem

.. math::
  \begin{align}
    \Delta u_1 &= f(\vec{x}) \\
    \frac{\partial u_1}{\partial \vec{n}} &= g(\vec{x}) \text{ on } \partial \omega
  \end{align}

within :math:`\omega`. :math:`u_2` is required to solve the Laplace equation while fixing the continuity condition on :math:`u`. This condition is met by a double-layer potential

.. math::
  \begin{align}
    u_2 = \int_{\partial \omega} u_1 \frac{\partial}{\partial \vec{n}} \frac{1}{|\vec{x} - \vec{x}'|} \dx
  \end{align}

References:
  [1] fredkin
