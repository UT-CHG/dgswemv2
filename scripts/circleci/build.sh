#!/bin/bash

cmake --version && g++ --version

cd $HOME/project
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=/home/circleci/install \
      -DSET_VERBOSE=ON -DUSE_OMPI=On -DUSE_HPX=On -DBUILD_EXAMPLES=On ..
make -j3
