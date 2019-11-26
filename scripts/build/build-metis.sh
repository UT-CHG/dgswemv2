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

METIS_BUILD_PATH="${BUILD_PATH}/metis"

if [ "$1" == "clean" ]; then
    clean_up $METIS_BUILD_PATH
fi

try_loading_modules $MODULES

if [ ! -d "$METIS_BUILD_PATH" ]; then
    set -e
    mkdir -p ${METIS_BUILD_PATH}
    cd ${METIS_BUILD_PATH}

    { #try pulling from metis website
        wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
    } || { #otherwise pull from the fedora mirror
        wget http://pkgs.fedoraproject.org/repo/pkgs/metis/metis-5.1.0.tar.gz/md5/5465e67079419a69e0116de24fce58fe/metis-5.1.0.tar.gz
    }

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
