#!/bin/bash

START_BOLD="tput bold"
END_BOLD="tput sgr0"

usage () {
    echo "usage: $0 [options]"
    echo "options:"
    echo "    -h"
    echo "        shows this message"
    echo "    -c <config_file> "
    echo "        specify a config file."
    echo "        default config file is config.txt"
    echo "    clean"
    echo "        removes build directory."
    echo "    no-make"
    echo "        runs CMake, does not run make."
    echo "    reinstall"
    echo "        if already built, reinstalls into install dir."
    exit 1
}

CONFIGFILE=""

if [ "$#" -gt 2 ]; then usage; fi

if [ "$#" == 2 ] && [ "$1" == "-c" ]; then
    CONFIGFILE=$2
fi

if [ "$1" == "-h" ]; then usage; fi
if [ "$1" == "--help" ]; then usage; fi

# Fix me: handle more options
if [ "$#" -gt 0 ]; then
    if [ "$1" != "clean" ] && [ "$1" != "no-make" ] && [ "$1" != "reinstall"]; then
	echo "invalid option: $1"
	usage
    fi
fi

SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
source $SCRIPTPATH/util.sh

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
                 -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH}"
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
