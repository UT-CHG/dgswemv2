#!/bin/bash

cmake --version && g++ --version
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j4
