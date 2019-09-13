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
source SCRIPTPATH/util.sh

if [ -z "$CONFIGFILE" ]; then
    CONFIGFILE=${SCRIPTPATH}/config.txt
fi

if [ -z "$CONFIGFILE" ]; then
    echo "CONFIGFILE variable not set, exiting!"
    exit 1
fi

load_config_file $CONFIGFILE

METIS_BUILD_PATH="${BUILD_PATH}/metis"

if [ "$1" == "clean" ]; then
    clean_up $METIS_BUILD_PATH
fi

try_loading_modules $MODULES

if [ ! -d "$METIS_BUILD_PATH" ]; then
    set -e
    mkdir -p ${METIS_BUILD_PATH}
    cd ${METIS_BUILD_PATH}
    wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
    tar xf metis-5.1.0.tar.gz
    cd metis-5.1.0
    sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/g' include/metis.h
    sed -i 's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/g' include/metis.h
    make config prefix=${INSTALL_PATH} CC=${C_COMPILER} CXX=${CXX_COMPILER}
    make
    make install
else
    $START_BOLD
    echo "directory exists! please either run:"
    $END_BOLD
    echo "$0 clean"
    $START_BOLD
    echo "or continue the build process:"
    $END_BOLD
    echo "cd ${METIS_BUILD_PATH}"
    echo "make"
    echo "make install"
    exit 0
fi
