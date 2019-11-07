#!/bin/bash

START_BOLD="tput bold"
END_BOLD="tput sgr0"

SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
source $SCRIPTPATH/util.sh

parse_args "$@"

if [ -z "$CONFIGFILE" ]; then
    CONFIGFILE=${SCRIPTPATH}/config.txt
fi

if [ -z "$CONFIGFILE" ]; then
    echo "CONFIGFILE variable not set, exiting!"
    exit 1
fi

load_config_file $CONFIGFILE

BLAZE_BUILD_PATH="$BUILD_PATH/blaze"

if [ "$1" == "clean" ]; then
    clean_up $BLAZE_BUILD_PATH
fi

try_loading_modules $MODULES

if [ ! -d "$BLAZE_BUILD_PATH" ]; then
    set -e
    mkdir -p $BLAZE_BUILD_PATH
    cd $BLAZE_BUILD_PATH
    git clone https://bitbucket.org/blaze-lib/blaze.git
    cd blaze
    BLAZE_REPO_PATH=$(pwd)
    mkdir build
    cd build
    CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
                 -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH}
                 -DBLAZE_SHARED_MEMORY_PARALLELIZATION=OFF"
    if [ -v CXX_COMPILER ]; then
    CMAKE_FLAGS="$CMAKE_FLAGS \
                 -DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    fi
    if [ -v C_COMPILER ]; then
    CMAKE_FLAGS="$CMAKE_FLAGS \
                 -DCMAKE_C_COMPILER=${C_COMPILER}"
    fi

    CMD="cmake ${CMAKE_FLAGS} ${BLAZE_REPO_PATH}"
    echo "CMD = $CMD"

    $CMD

    make -k -j$NUM_BUILDCORES
    make install
else
    $START_BOLD
    echo "directory exists! please either run:"
    $END_BOLD
    echo "$0 clean"
    $START_BOLD
    echo "or continue the build process:"
    $END_BOLD
    echo "cd ${BLAZE_BUILD_PATH}"
    echo "make"
    echo "make install"
    exit 0
fi
