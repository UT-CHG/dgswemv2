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
JEMALLOC_BUILD="${BUILD_PATH}/jemalloc"
if [ "$1" = "clean" ]; then
    clean_up $JEMALLOC_BUILD
fi

if [ ! -d $INSTALL_PATH ]; then
    echo "Creating install path..."
    mkdir -p ${INSTALL_PATH}
fi

try_loading_modules $MODULES

if [ ! -d ${JEMALLOC_BUILD} ]; then
    set -e
    mkdir -p ${JEMALLOC_BUILD}
    cd ${JEMALLOC_BUILD}
    wget https://github.com/jemalloc/jemalloc/releases/download/3.6.0/jemalloc-3.6.0.tar.bz2
    tar xf jemalloc-3.6.0.tar.bz2
    cd jemalloc-3.6.0
    ./configure --prefix=$INSTALL_PATH CC=${C_COMPILER} CXX=${CXX_COMPILER}
    make -j install
    exit
else
    set -e
    cd ${JEMALLOC_BUILD}/jemalloc-3.6.0
    make install
    exit
fi
