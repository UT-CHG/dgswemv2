#!/bin/bash

cmake --version && g++ --version

cd $HOME/dgswemv2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=/home/ubuntu/install ..
make -j4
