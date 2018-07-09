#!/bin/bash

cmake --version && g++-5 --version

#clone yaml-cpp
cd /usr
git clone https://github.com/jbeder/yaml-cpp.git
cd /usr/dgswemv2/scripts/
#Use an implicit newline in echo to simulate hitting the enter key
echo | build/build-yaml-cpp.sh -c docker/docker.config.txt

echo | build/build-boost.sh -c docker/docker.config.txt

echo | build/build-jemalloc.sh -c docker/docker.config.txt

echo | build/build-metis.sh -c docker/docker.config.txt

#Build linear algebra libraries
echo | build/build-eigen.sh -c docker/docker.config.txt

echo | build/build-blaze.sh -c docker/docker.config.txt

#hpx
cd /usr
git clone https://github.com/STEllAR-GROUP/hpx.git
cd /usr/dgswemv2/scripts/
echo | build/build-hpx.sh -c docker/docker.config.txt

#clean up
cd /usr
rm -rf yaml-cpp
rm -rf hpx