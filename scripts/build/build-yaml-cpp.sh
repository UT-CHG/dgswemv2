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

if [ -d "${YAML_CPP_REPO_PATH}" ]; then
    echo "found yaml-cpp repo path: ${YAML_CPP_REPO_PATH}"
else
    echo "could not find yaml-cpp repo path: ${YAML_CPP_REPO_PATH}"
    echo "clone yaml-cpp to ${YAML_CPP_REPO_PATH}?"
    echo "Press enter to continue with these settings (or ctrl-c to exit)"
    read answer

    git clone https://github.com/jbeder/yaml-cpp.git ${YAML_CPP_REPO_PATH}
fi

if [ "$1" == "clean" ]; then
    clean_up $YAML_CPP_BUILD_PATH
fi

if [ "$1" == "reinstall" ]; then
    if [ -d ${YAML_CPP_BUILD_PATH} ]; then
	cd ${YAML_CPP_BUILD_PATH}
	make install
	if [ "$#" != 0 ]; then
	    $START_BOLD; echo "make install may have failed!"; $END_BOLD;
	fi
    fi
fi

echo "Build type is: ${BUILD_TYPE}"
echo "Build path is: ${YAML_CPP_BUILD_PATH}"

try_loading_modules $MODULES

if [ ! -d "$YAML_CPP_BUILD_PATH" ]; then
    echo "Creating build path..."
    mkdir -p "$YAML_CPP_BUILD_PATH"
    cd "$YAML_CPP_BUILD_PATH"

    ############################################
    # MODIFY CMAKE FLAGS AS YOU WISH
    #

    CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
                 -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH}\
                 -DYAML_CPP_BUILD_TESTS=OFF"
    if [ -v CXX_COMPILER ]; then
    CMAKE_FLAGS="$CMAKE_FLAGS \
                 -DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    fi
    if [ -v C_COMPILER ]; then
    CMAKE_FLAGS="$CMAKE_FLAGS \
                 -DCMAKE_C_COMPILER=${C_COMPILER}"
    fi

    CMD="cmake ${CMAKE_FLAGS} $YAML_CPP_REPO_PATH"
    echo "CMD = $CMD"

    $CMD

    #
    ############################################


    # cleanup on failure
    rc=$?
    if [[ $rc != 0 ]] ; then
	echo "build failure, deleting build path"
        cd "$SCRIPTPATH"
        rm -rf "$YAML_CPP_BUILD_PATH"
        exit $rc
    fi

else
    $START_BOLD
    echo "directory exists! please either run:"
    $END_BOLD
    echo "$0 clean"
    $START_BOLD
    echo "or continue the build process:"
    $END_BOLD
    echo "cd ${YAML_CPP_BUILD_PATH}"
    echo "nice make -k -j${NUM_BUILDCORES} || exit 1"
    echo "make install"
    exit 0
fi

# the actual build command
if [ "$1" != "no-make" ]; then
    cd "$YAML_CPP_BUILD_PATH"
    MAKECMD="make -k -j$NUM_BUILDCORES"
    INSTALLCMD="make install"

    $MAKECMD
    $INSTALLCMD
fi
