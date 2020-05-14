#!/bin/bash

cmake --version && g++ --version

#clone yaml-cpp
cd /usr
git clone https://github.com/jbeder/yaml-cpp.git
cd yaml-cpp
git checkout tags/yaml-cpp-0.6.3

cd /usr/dgswemv2/scripts/
#Use an implicit newline in echo to simulate hitting the enter key
echo | build/build-yaml-cpp.sh -c docker/docker.config.txt

echo | build/build-metis.sh -c docker/docker.config.txt

#Build linear algebra libraries
echo | build/build-eigen.sh -c docker/docker.config.txt

echo | build/build-blaze.sh -c docker/docker.config.txt

#clean up
cd /usr
rm -rf yaml-cpp