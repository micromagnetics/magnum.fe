Getting Started
===============

As a first example the standard problem #4 [MuMag4]_, proposed by the MuMag group is computed with magnum.fe. Since magnum.fe is a Python library, a simulation script is a Python script that imports magnum.fe. Thus, every magnum.fe simulation script starts with

.. code:: python

  from magnumfe import *

In the next step a mesh for the simulation is created. For simple geometries, the builtin meshing tools of dolfin and magnum.fe are a good choice. More complicated geometries can be meshed with external tools. magnum.fe supports a variety of mesh file formats through the Gmsh library. Depending on the method for the solution of the open-boundary demagnetization field problem the mesh is also required to include a cuboid shell, see :ref:`open-boundary`.

Here we use the hybrid FEM-BEM method, see :ref:`fem-bem`, for the computation of the demagnetization field. Hence a simple cuboid mesh as provided by the dolfin meshing tools is sufficient for the solution of the standard problem #4, that requires a cuboid of size :math:`500 \times 125 \times 3` nm. The mesh is constructed symmetrically around the coordinate origin and scaled by :math:`10^9`.

.. code:: python

  mesh = BoxMesh(-500.0/2, -125.0/2, -3.0/2, 500.0/2, 125.0/2, 3.0/2, 100, 25, 1)

A magnetization configuration that is known to relax quickly in a magnetic s-state is defined by an analytical expression.

.. code:: python

  m_start = Expression(("cos(alpha)", "sin(alpha)", "0.0"), alpha=Expression("fabs(pi*x[0]/1e3)"))

A simulation state is created and initialized with material paramters and the start magnetization m_start. Note the the :code:`scale` parameter is set to :math:`10^{-9}` since the problem size is given in nanometers. This scaling is advised to avoid numerical artifacts.

.. code:: python

  state = State(mesh, material = Material.py(), scale = 1e-9, m = m_start)

For the relaxation of the system, an LLG-solver object is created that includes the exchange field and demagnetization field as only effective field contributions.

.. code:: python

  llg = LLGAlougesProject([ExchangeField(), DemagField("FK")])

The system is relaxed by setting the damping of the material to :math:`\alpha = 1` and performing a number of integration steps.

.. code:: python

  state.material.alpha = 1.0
  for i in range(200): llg.step(state, 2e-11)

In order to switch the magnetization as required by the standard problem #4, the damping is reduced to :math:`\alpha = 0.02`. A new solver object is created that includes an external-field.

.. code:: python

  state.material.alpha = 0.02

  llg = LLGAlougesProject([
      ExternalField(Constant((-24.6e-3/Constants.mu0, +4.3e-3/Constants.mu0, 0.0))),
      ExchangeField(),
      DemagField("FK")
  ])

The time loop for the solution of the LLG has to be programmed explicitly by now. Also the logging of the averaged magnetization is realized directly in Python. Note that an LLG step can either be performed by calling :code:`step` on the :class:`llg` object as is done for the relaxation process, or by calling :code:`step` on the :class:`state` object. In contrast to the first method, the latter method increases the time variable :class:`t` of the state.

.. code:: python

  # open logfile
  logfile = open("sp4.dat", "w", 0)

  # initialize time variables
  dt, T = 2e-13, 1e-9

  # loop, loop, loop
  for i in range(int(T / dt)):
    
    # write scalar information
    logfile.write("%f %f %f %f\n" % ((state.t*1e9,) + state.m.average()))

    # calculate next step
    state.step(llg, dt)

  logfile.close()

Complete code
+++++++++++++

.. code:: python

  from magnumfe import *

  #######################################
  #### GENERATE MESH WITH SHELL
  #######################################

  mesh = BoxMesh(-500.0/2, -125.0/2, -3.0/2, 500.0/2, 125.0/2, 3.0/2, 100, 25, 1)

  #######################################
  #### RELAX SYSTEM TO S-STATE
  #######################################

  # define start magnetization
  m_start = Expression(("cos(alpha)", "sin(alpha)", "0.0"), alpha=Expression("fabs(pi*x[0]/1e3)"))

  state   = State(mesh, material = Material.py(), scale = 1e-9, m = m_start)
  llg     = LLGAlougesProject([ExchangeField(), DemagField("FK")])

  state.material.alpha = 1.0
  for i in range(200): llg.step(state, 2e-11)

  #######################################
  #### SIMULATE SWITCHING
  #######################################

  state.material.alpha = 0.02

  llg = LLGAlougesProject([
      ExternalField(Constant((-24.6e-3/Constants.mu0, +4.3e-3/Constants.mu0, 0.0))),
      ExchangeField(),
      DemagField("FK")
  ])

  logfile = open("sp4.dat", "w", 0)
  dt, T = 2e-13, 1e-9

  for i in range(int(T / dt)):
    
    # write scalar information
    logfile.write("%.10f %f %f %f\n" % ((state.t*1e9,) + state.m.average()))

    # calculate next step
    state.step(llg, dt)

  logfile.close()
 
Run the Simulation
++++++++++++++++++

Since the simulation file is a simple Python script it is run with the Python interpreter. Save the above program to a file called `sp4.py` and run

.. code::

  $ python sp4.py

on the command line.

More Examples
+++++++++++++

More examples can be found in the :code:`examples` directory of the magnum.fe source tree.


.. [MuMag4] µMAG Standard Problem #4, http://www.ctcms.nist.gov/~rdm/std4/spec4.html

