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
    exit 1
}

CONFIGFILE=""

if [ "$#" -gt 2 ]; then usage; fi

if [ "$#" == 2 ] && [ "$1" == "-c" ]; then
    CONFIGFILE=$2
fi

if [ "$1" == "-h" ]; then usage; fi
if [ "$1" == "--help" ]; then usage; fi

if [ "$#" -gt 0 ]; then
    if [ "$1" != "clean" ] && [ "$1" != "no-make" ]; then
	echo "invalid option: $1"
	usage
    fi
fi

SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
export $SCRIPTPATH/util.sh

if [ -z "$CONFIGFILE" ]; then
    CONFIGFILE=${SCRIPTPATH}/config.txt
fi

if [ -z "$CONFIGFILE" ]; then
    echo "CONFIGFILE variable not set, exiting!"
    exit 1
fi

load_config_file $CONFIGFILE

# Done setting up variables.

BOOST_BUILD="${BUILD_PATH}/boost"
if [ "$1" = "clean" ]; then
    clean_up $BOOST_BUILD
fi

if [ ! -d $INSTALL_PATH ]; then
    echo "Creating install path..."
    mkdir -p ${INSTALL_PATH}
fi

try_loading_modules $MODULES

if [ ! -d ${BOOST_BUILD} ]; then
    set -e
    mkdir -p ${BOOST_BUILD}
    cd ${BOOST_BUILD}
    wget 'http://downloads.sourceforge.net/project/boost/boost/1.63.0/boost_1_63_0.tar.bz2'
    tar xf boost_1_63_0.tar.bz2
    cd boost_1_63_0
    ./bootstrap.sh --prefix="$INSTALL_PATH"
    ./b2 -j${NUM_BUILDCORES} install
    exit
else
    set -e
    cd ${BOOST_BUILD}/
    make install
    exit
fi
