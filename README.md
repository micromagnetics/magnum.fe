magnum.fe
=========

Installation
------------

### Prerequisites
magnum.fe requires the following software/libraries to be installed:

* FEniCS >= 1.4
* gmsh >= 2.8.0 (headers and library)
* dev-version of cbc.block (optional)
* BEM++ 2.0.1 (optional)
* cmake >= 2.8

#### Install dependencies in Ubuntu 14.04
Install FEniCS

    $ sudo add-apt-repository ppa:fenics-packages/fenics
    $ sudo apt-get update
    $ sudo apt-get install fenics

Install Gmsh with headers

    $ sudo apt-get install libgmsh-dev

Install CBC.Block

    $ bzr branch lp:cbc.block
    $ cd cbc.block
    $ sudo python setup.py install

For installation of BEM++, see http://www.bempp.org/. magnum.fe requires the H-matrix library AHMED 1.0 that is optional for BEM++, see http://bebendorf.ins.uni-bonn.de/AHMED.html.

### Build and install
To build magenum.fe with cmake do

    cd /path/to/magnum.fe
    mkdir build
    cd build
    cmake ..
    sudo make install

### Create docker container
magnum.fe can be virtualized with docker (http://www.docker.com) for easy deployment to other machines.
In order to create a docker container run

    $ cd /path/to/magnum.fe/docker
    $ cp /path/to/AHMED-1.0.tar.gz .
    $ sudo docker build .

Note that you have to place the AHMED 1.0 tarball in the `docker` directory for a successfull build.

### Test
You can test your installation by running the unit tests

    cd tests
    python run_all.py

License
-------
Released under the GNU General Public License 3 (included).
