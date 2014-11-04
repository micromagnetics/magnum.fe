.. only:: html

  .. image:: /images/logo.svg

Welcome to magnum.fe's documentation!
=====================================

magnum.fe is a finite-element package for the solution of dynamical micromagnetic problems. It is written in Python and C++ and largely relies on the multi purpose finite-element library `FEniCS <http://fenicsproject.org/>`_. For the solution of open-boundary problems a hybrid FEM/BEM method is implemented that uses the `BEM++ <http://www.bempp.org/>`_ library. magnum.fe is free software and can be extended to your needs.

Features
++++++++
* Integration of the Landau-Lifshitz-Gilbert Equation
* Demagnetization-field computation via shell transformation or hybrid FEM/BEM method
* Oersted-field computation via shell transformation or hybrid FEM/BEM method
* Solution of the spin-diffusion model

Download
++++++++
Visit our `GitHub page <https://github.com/micromagnetics/magnum.fe>`_ to download magnum.fe.

Citing
++++++
If you use results from magnum.fe for scientific publications, we would be very grateful if you would cite the following paper:

  Abert, C., Exl, L., Bruckner, F., Drews, A., & Suess, D. (2013). magnum. fe: A micromagnetic finite-element simulation code based on FEniCS. Journal of Magnetism and Magnetic Materials, 345, 29-35.

License and Disclaimer
++++++++++++++++++++++
magnum.fe is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

magnum.fe is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

Alternatives
++++++++++++
You might also want to have a look at these open-source finite-element micromagnetic codes:

* Nmag (http://nmag.soton.ac.uk/nmag/)
* magpar (http://www.magpar.net/)

You code is not listed here, but should be? Drop us a line (claas.abert@tuwien.ac.at).

.. toctree::
   :maxdepth: 2

   installation
   getting_started
   state
   meshing
   llg
   spin_diffusion
   open_boundary
   field_terms

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
