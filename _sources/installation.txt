Installation
============

Prerequisites
+++++++++++++
magnum.fe requires the following software/libraries to be installed:

* FEniCS >= 1.5
* CMake >= 2.8
* SWIG >= 2.0
* G++ >= 4.0
* dev-version of CBC.Block (optional)
* BEM++ 2.0.3 python-pseudoinverse branch (optional)

Install from source
+++++++++++++++++++

The current version of magnum.fe can be downloaded from GitHub via

.. code::

  $ git clone https://github.com/micromagnetics/magnum.fe.git

To build und install magnum.fe with CMake do

.. code::

  $ cd /path/to/magnum.fe
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
  $ sudo make install

Install dependencies in Ubuntu 14.04
++++++++++++++++++++++++++++++++++++
Install FEniCS

.. code::

  $ sudo add-apt-repository ppa:fenics-packages/fenics
  $ sudo apt-get update
  $ sudo apt-get install fenics

Install CBC.Block

.. code::

  $ git clone https://bitbucket.org/fenics-apps/cbc.block.git
  $ cd cbc.block
  $ sudo python setup.py install

For installation of BEM++, see http://www.bempp.org/.

Create docker container
+++++++++++++++++++++++
magnum.fe can be virtualized with docker (http://www.docker.com) for easy deployment to other machines.
In order to create a docker container run

.. code::

  $ cd /path/to/magnum.fe
  $ sudo docker build .

This procedure should download all necessary software and install magnum.fe into the container.
You need Docker >= 1.0 to build the container.
