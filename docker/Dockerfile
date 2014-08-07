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
    git

RUN cd /tmp && \
    git clone -b dev https://github.com/c-abird/magnum.fe.git && \
    cd magnum.fe && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

# Cleanup to save space
RUN apt-get clean && \ 
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the HOME environment variable, otherwise import dolfin crashes
RUN echo "/root" > /etc/container_environment/HOME

# Set LIBGL_ALWAYS_INDIRECT to suppress libGL warning message.
RUN echo "y" > /etc/container_environment/LIBGL_ALWAYS_INDIRECT

CMD ["/sbin/my_init"]
