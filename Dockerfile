# Builds a Docker image with the official FEniCS PPA packages 
# and magnum.fe.
#
# Authors: 
# Claas Abert <claas.abert@tuwien.ac.at>
#
# Based on the work of
# Lizao Li <lixx1445@umn.edu>
# Jack S. Hale <jack.hale@uni.lu>
# (https://bitbucket.org/garth-wells/fenics-virtual/src/4efc1c04b8c83b50497d2c1e7d239019bd429569/docker/stable-ppa/Dockerfile)
FROM phusion/baseimage:0.9.12
#FROM fenicsproject/stable
MAINTAINER Claas Abert <claas.abert@tuwien.ac.at>

ENV HOME /root

# Install add-apt-repository
RUN apt-get -qq update && \
    apt-get -qqy install python-software-properties

# Install the basic environment and fenics from the PPA
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qqy install xauth fenics ipython

# Install basic tools and compilers
RUN apt-get -qq update && \
    apt-get -qqy install cmake swig g++ git gfortran
    
# Install CBC.Block
RUN cd /tmp && \
    git clone https://bitbucket.org/micromagnetics/cbc.block.git && \
    cd cbc.block && \
    python setup.py install

# Install Bempp
RUN apt-get -qq update && \
    apt-get -qqy install python-m2crypto

RUN cd /tmp && \
    git clone -b features/python-pseudoinverse https://github.com/micromagnetics/bempp.git

COPY docker/bempp_setup.cfg /tmp/bempp/

RUN cd /tmp/bempp && \
    python bempp_setup.py -b bempp_setup.cfg && \
    python bempp_setup.py -c bempp_setup.cfg && \
    python bempp_setup.py -i all bempp_setup.cfg 

ENV PYTHONPATH $HOME/bempp/python:$PYTHON_PATH



# Install GMSH
RUN apt-get -qq update && \
    apt-get -qqy install libgmsh-dev

# Install magnum.fe
COPY . /tmp/magnum.fe

RUN cd /tmp/magnum.fe && \
    rm -rf build && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    ldconfig

# Remove magnum.fe source files
RUN cd /usr/local/lib/python2.7/dist-packages/magnumfe && \
    find . -name "*.py"| xargs rm

# Cleanup to save space
RUN apt-get clean && \ 
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the HOME environment variable, otherwise import dolfin crashes
RUN echo "/root" > /etc/container_environment/HOME

# Set LIBGL_ALWAYS_INDIRECT to suppress libGL warning message.
RUN echo "y" > /etc/container_environment/LIBGL_ALWAYS_INDIRECT

CMD ["/sbin/my_init"]
