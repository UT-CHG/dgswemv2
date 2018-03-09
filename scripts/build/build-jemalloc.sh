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

if ["$#" -gt 0]; then
    if [ "$1" != "clean" ] && [ "$1" != "no-make" ]; then
	echo "invalid option: $1"
	usage
    fi
fi

SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
if [ -z "$CONFIGFILE" ]; then
    CONFIGFILE=${SCRIPTPATH}/config.txt
fi

if [ -z "$CONFIGFILE" ]; then
    echo "CONFIGFILE variable not set, exiting!"
    exit 1
fi

echo -n "config file: "
$START_BOLD
echo "$CONFIGFILE"
$END_BOLD
while read line; do

    if [[ "$line" =~ "#" ]]; then
	echo -e "    \e[31m${line}\e[0m"
    else    
	VAR=`echo $line | awk -F= '{print $1}'`
	VALUE=`echo $line | awk -F= '{print $2}'`
	if [ -z "$VAR" ] || [ -z "$VALUE" ]; then
	    echo "    $line"
	else
	    echo -n "    "
	    echo -e "\e[33m${VAR}\e[0m=\e[37m${VALUE}\e[0m"
	fi
    fi
done < $CONFIGFILE

echo "Press enter to continue with these settings (or ctrl-c to exit)."
read answer

if [ -s "$CONFIGFILE" ]; then
    source $CONFIGFILE
else
    echo "config file $CONFIGFILE either empty or non-existent"
    exit 1
fi

if [ -z "$VERBOSE" ]; then
    VERBOSE="false"
fi

if [ "$VERBOSE" = "true" ]; then
    echo "verbose mode"
    set -x
else
    set +x
fi

# Done setting up variables.
JEMALLOC_BUILD="${BUILD_PATH}/jemalloc"
if [ "$1" = "clean" ]; then
    CLEAN_CMD="rm -rf ${JEMALLOC_BUILD}"

    $START_BOLD
    echo "$0 clean:"
    $END_BOLD
    echo "about to execute:"
    echo "${CLEAN_CMD}"
    $START_BOLD
    echo "is this okay? [y/N]"
    $END_BOLD
    read answer
    if echo "$answer" | grep -iq "^y" ;then
	echo "removing build directory ${JEMALLOC_BUILD}"
	$CLEAN_CMD
    else
	echo "doing nothing."
    fi
    exit 0
fi

if [ ! -d $INSTALL_PATH ]; then
    echo "Creating install path..."
    mkdir -p ${INSTALL_PATH}
fi

for module in $MODULES; do
    module load $module
done

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
