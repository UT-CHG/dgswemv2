#!/bin/bash

cmake --version && g++ --version

#clone yaml-cpp
cd $HOME
git clone https://github.com/jbeder/yaml-cpp.git

cd $HOME/dgswemv2/scripts/
#Use an implicit newline in echo to simulate hitting the enter key
echo | build/build-yaml-cpp.sh -c circleci/test.config.txt

cd $HOME/dgswemv2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
