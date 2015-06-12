.. module:: magnumfe

Meshing
=======

magnum.msh
++++++++++

.. figure:: /images/mesh.png
   :width: 100%
   :align: center
   
   Mesh of a sphere including cuboid shell for the solution of open-boundary problems.
   (a) Explosion view
   (b) Cross section

magnum.fe requires mesh files in the FEniCS XML format. Domain information may either be included in the mesh file or defined in additional XML files as :class:`dolfin.MeshFuntion`. To convert your existing mesh files, you can either use the FEniCS tool :code:`dolfin-convert` or the magnum.fe tool `magnum.msh <https://github.com/micromagnetics/magnum.msh>`_.

Application of the shell-transformation method for the solution of open-boundary problems requires meshes with properly defined shell elements, see figure. Given a mesh of the magnetic sample, magnum.msh is able to automatically create a mesh including the required shell elements.

:class:`WrappedMesh`
++++++++++++++++++++

.. autoclass:: WrappedMesh
   :members:
