Installation
============

Prerequisites
-------------
magnum.fe requires the following software/libraries to be installed:

* FEniCS >= 1.4
* Gmsh >= 2.8.0 (headers and library)
* CMake >= 2.8
* SWIG >= 2.0
* G++ >= 4.0
* dev-version of CBC.Block (optional)
* BEM++ 2.0.1 (optional)

Install from source
-------------------

The current version of magnum.fe can be downloaded from GitHub via

.. code::

  $ git clone https://github.com/micromagnetics/magnum.fe.git

To build und install magenum.fe with CMake do

.. code::

  $ cd /path/to/magnum.fe
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
  $ sudo make install

Install dependencies in Ubuntu 14.04
------------------------------------
Install FEniCS

.. code::

  $ sudo add-apt-repository ppa:fenics-packages/fenics
  $ sudo apt-get update
  $ sudo apt-get install fenics

Install Gmsh with headers and build dependencies

.. code::

  $ sudo apt-get install libgmsh-dev g++ swig cmake

Install CBC.Block

.. code::

  $ git clone https://bitbucket.org/fenics-apps/cbc.block.git
  $ cd cbc.block
  $ sudo python setup.py install

For installation of BEM++, see http://www.bempp.org/. magnum.fe requires the H-matrix library AHMED 1.0 that is optional for BEM++, see http://bebendorf.ins.uni-bonn.de/AHMED.html.

Create docker container
-----------------------
magnum.fe can be virtualized with docker (http://www.docker.com) for easy deployment to other machines.
In order to create a docker container run

.. code::

  $ cd /path/to/magnum.fe/docker
  $ cp /path/to/AHMED-1.0.tar.gz .
  $ sudo docker build .

You need Docker >= 1.0 to build the container. Note that you have to place the AHMED 1.0 tarball in the `docker` directory for a successfull build (http://bebendorf.ins.uni-bonn.de/AHMED.html).
