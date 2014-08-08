# Builds a Docker image with the official FEniCS PPA packages 
# and magnum.fe.
#
# Authors: 
# Claas Abert <claas.abert@tuwien.ac.at>
FROM phusion/baseimage:0.9.12
MAINTAINER Claas Abert <claas.abert@tuwien.ac.at>

ENV HOME /root

# Install add-apt-repository
RUN apt-get -qq update && \
    apt-get -qqy install python-software-properties

# Install the basic environment and fenics from the PPA
RUN add-apt-repository -y ppa:fenics-packages/fenics && \
    apt-get -qq update && \
    apt-get -qqy install xauth fenics ipython

# Install Compiler and magnum.fe dependencies 
RUN apt-get update && \
    apt-get install -yy \
    cmake \
    swig \
    g++ \
    libgmsh-dev \
    git \
    gfortran
    
# Install CBC.Block
#RUN apt-get update && \
#    apt-get install -yy bzr
#RUN cd /tmp && \
#    bzr branch lp:cbc.block && \
#    cd cbc.block && \
#    python setup.py install

# Install Bempp
RUN cd /tmp && \
    git clone -b release_2.0 https://github.com/bempp/bempp.git

COPY bempp_setup.cfg /tmp/bempp/bempp_setup.cfg
COPY AHMED-1.0.tar.gz /tmp/AHMED-1.0.tar.gz

RUN cd /tmp/bempp && \
    python bempp_setup.py -b bempp_setup.cfg && \
    python bempp_setup.py -c bempp_setup.cfg && \
    python bempp_setup.py -i all bempp_setup.cfg 

ENV PYTHONPATH $HOME/bempp/python:$PYTHON_PATH

# Compile and install magnumfe
RUN cd /tmp && \
    git clone -b dev https://github.com/c-abird/magnum.fe.git && \
    cd magnum.fe && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    ldconfig

# Cleanup to save space
RUN apt-get clean && \ 
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the HOME environment variable, otherwise import dolfin crashes
RUN echo "/root" > /etc/container_environment/HOME

# Set LIBGL_ALWAYS_INDIRECT to suppress libGL warning message.
RUN echo "y" > /etc/container_environment/LIBGL_ALWAYS_INDIRECT

CMD ["/sbin/my_init"]
