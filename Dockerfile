FROM ubuntu:18.04

################################################################################
# Docker basics
#
# To build the container, run for <DGSWEMV2_ROOT>:
#
# docker build .
#
# You can use the `--no-cache` argument to build the container from scratch
# and `-t <tag> to build the image with a specific tag
#
# To push the image to dockerhub, run:
#
# docker push <dockerhub_id>/dgswemv2:<tag>
#
# Finally to run the container locally, run:
#
#
# docker run -it <dockerhub_id>/dgswemv2:<tag>
#
# to start an interactive shell with the given docker image

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
    apt-get -y install python-dev                              && \
    apt-get -y install python-pip                              && \
    apt-get -y install git                                     && \
    pip install numpy

# Set the default working directory to dgswemv2 dir
WORKDIR /usr/dgswemv2

# Load source into container
COPY . /usr/dgswemv2

# Setup project space
RUN /usr/dgswemv2/scripts/docker/set-up-dependencies.sh
#    cd dgswemv2                                            && \

# Final project setup
#RUN sh setup.sh

CMD /bin/bash