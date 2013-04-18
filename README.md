magnum.fe
=========

Installation
------------

### Prerequisites
magnum.fe requires the following software/libraries to be installed:

* dev-version of FEniCS
* dev-version of cbc.block
* gmsh (headers and library)
* cmake

magnum.fe is tested with dolfin build 7614, but should also work with more recent versions.

#### Install dev-version of dolfin
The dev-version of FEniCS is best installed with dorsal (https://launchpad.net/dorsal).

    bzr branch lp:dorsal
    cd dorsal

Set `STABLE_BUILD=false` in the dorsal.cfg. Then run
    
    ./dorsal.sh

and follow the instructions.

#### Install dev-version of cbc.block
The dev-version of cbc.block can be installed by checking out the current version from the launchpad repository

    bzr branch lp:cbc.block
    cd cbc.block
    sudo python setup.py install

### Build
To build magenum.fe with cmake do

    cd /path/to/magnum.fe
    mkdir build
    cd build
    cmake ..
    make

Right now there is no install rule for make. There are however symlinks in the `site-packages/magnumfe` directory which point to the compiled swig wrapper code in the `build` directory. Thus you can use magnum.fe by adding the `site-packages` directory to your `PYTHONPATH`

    export PYTHONPATH=/path/to/magnum.fe/site-packages

### Test
You can test your installation by running the unit tests

    cd tests
    python test.py

License
-------
Released under the GNU General Public License 3 (included).
