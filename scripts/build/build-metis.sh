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
if ["$#" -gt 0]; then
    if [ "$1" != "clean" ] && [ "$1" != "no-make" ] && [ "$1" != "reinstall"]; then
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

if [ "$1" == "clean" ]; then
    CLEAN_CMD="rm -rf $METIS_BUILD_PATH"

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
	echo "removing build directory ${METIS_BUILD_PATH}"
	$CLEAN_CMD
    else
	echo "doing nothing."
    fi
    exit 0
fi

echo "MODULES = $MODULES"

module purge

for module in $MODULES; do
    module load $module
done

module list

#echo "git branch:"
#git show-branch >> $LOGFILE
#echo "git hash:"
#git rev-parse HEAD >> $LOGFILE
#echo "git status:"
#git status >> $LOGFILE

#echo "*** config file: ***"
#cat $CONFIGFILE >> $LOGFILE
#echo "*** end config file ***"

#echo "MODULES = $MODULES" >> $LOGFILE

METIS_BUILD="${BUILD_PATH}/metis"
if [ ! -d "$METIS_BUILD_PATH" ]; then
    set -e
    mkdir -p ${METIS_BUILD}
    cd ${METIS_BUILD}
    wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
    tar xf metis-5.1.0.tar.gz
    cd metis-5.1.0
    sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/g' include/metis.h
    sed -i 's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/g' include/metis.h
    make config prefix=${INSTALL_PATH}
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
