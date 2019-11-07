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
    wget 'https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2'
    tar xf boost_1_69_0.tar.bz2
    cd boost_1_69_0
    ./bootstrap.sh --prefix="$INSTALL_PATH"
    ./b2 -j${NUM_BUILDCORES} install
    exit
else
    set -e
    cd ${BOOST_BUILD}/
    make install
    exit
fi
