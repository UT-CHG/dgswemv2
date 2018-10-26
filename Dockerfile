FROM ubuntu:18.04

# github.com/bremerm31
MAINTAINER Max Bremer

# Install all tools + libraries needed for the project
RUN apt-get update                                             && \
    apt-get -y install apt-utils                               && \
    apt-get -y install ssh                                     && \
    apt-get -y install build-essential                         && \
    apt-get -y install cmake                                   && \
    apt-get -y install libz-dev                                && \
    apt-get -y install openmpi-bin libopenmpi-dev              && \
    apt-get -y install gfortran                                && \
    apt-get -y install libblas-dev liblapack-dev               && \
    apt-get -y install libhwloc-dev                            && \
    apt-get -y install libboost-all-dev                        && \
    apt-get -y install python-dev                              && \
    apt-get -y install python-pip                              && \
    apt-get -y install git                                     && \
    pip install numpy

# Set the default working directory to dgswemv2 dir
WORKDIR /usr/dgswemv2

# Load source into container
COPY . /usr/dgswemv2

# Setup project space
#RUN /usr/dgswemv2/scripts/docker/set-up-dependencies.sh
#    cd dgswemv2                                            && \


# Final project setup
#RUN sh setup.sh

CMD /bin/bash