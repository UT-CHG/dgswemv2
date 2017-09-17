#!/bin/bash

cmake --version && g++ --version

#clone yaml-cpp
cd $HOME
git clone https://github.com/jbeder/yaml-cpp.git

cd $HOME/dgswemv2/scripts/
#Use an implicit newline in echo to simulate hitting the enter key
echo | build/build-yaml-cpp.sh -c circleci/test.config.txt

echo | build/build-boost.sh -c circleci/test.config.txt

echo | build/build-jemalloc.sh -c circleci/test.config.txt

echo | build/build-metis.sh -c circleci/test.config.txt

#hpx
cd $HOME
git clone https://github.com/STEllAR-GROUP/hpx.git
cd $HOME/dgswemv2/scripts/
echo | build/build-hpx.sh -c circleci/test.config.txt
