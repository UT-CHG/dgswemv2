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

if [ "$1" == "reinstall" ]; then
    $START_BOLD
    echo "WARNING! reinstalling old files if you've changed the config file might be a bad idea"
fi

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

if [ -d "${HPX_REPO_PATH}" ]; then
    echo "found HPX repo path: ${HPX_REPO_PATH}"
else
    echo "could not find HPX repo path: ${HPX_REPO_PATH}!"
    exit 1
fi

if [ "$1" == "clean" ]; then
    CLEAN_CMD="rm -rf $HPX_BUILD_PATH"

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
	echo "removing build directory ${HPX_BUILD_PATH}"
	$CLEAN_CMD
    else
	echo "doing nothing."
    fi
    exit 0
fi

if [ "$1" == "reinstall" ]; then
    if [ -d ${HPX_BUILD_PATH} ]; then
	cd ${HPX_BUILD_PATH}	
	make install
	if [ "$#" != 0 ]; then
	    $START_BOLD; echo "make install may have failed!"; $END_BOLD;
	fi
    fi
fi

echo "Build type is: ${BUILD_TYPE}"
echo "Build path is: ${HPX_BUILD_PATH}"

# Add VTune module
if [ "$VTUNE" = "true" ]; then
    MODULES="${MODULES} ${VTUNE_MODULE}"
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

if [ ! -d "$HPX_BUILD_PATH" ]; then
    echo "Creating build path..."
    mkdir -p "$HPX_BUILD_PATH"
    cd "$HPX_BUILD_PATH"

    ############################################
    # MODIFY CMAKE FLAGS AS YOU WISH
    #

    MALLOC="jemalloc"
    IDLE_RATES="true"

    CMAKE_FLAGS="-DCMAKE_TOOLCHAIN_FILE=${HPX_REPO_PATH}/cmake/toolchains/Stampede2-gcc.cmake \
                 -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
                 -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} \
                 -DHPX_WITH_PARCELPORT_MPI=true \
                 -DHPX_WITH_MALLOC=$MALLOC \
                 -DHPX_WITH_THREAD_IDLE_RATES=${IDLE_RATES} \
                 -DHPX_WITH_CXX14=On"
    if [ $VTUNE = "true" ]; then
	CMAKE_FLAGS="${CMAKE_FLAGS} \
                 -DHPX_WITH_ITTNOTIFY=On \
                 -DAMPLIFIER_ROOT=${VTUNE_DIR}"
    fi

    CMD="cmake ${CMAKE_FLAGS} $HPX_REPO_PATH"
    echo "CMD = $CMD"

    $CMD

    #
    ############################################


    # cleanup on failure
    rc=$?
    if [[ $rc != 0 ]] ; then
	echo "build failure, deleting build path"
        cd "$SCRIPTPATH"
        rm -rf "$HPX_BUILD_PATH"
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
    echo "cd ${HPX_BUILD_PATH}"
    echo "nice make -k -j${NUM_BUILDCORES} || exit 1"
    echo "make install"
    exit 0
fi

# the actual build command
if [ "$1" != "no-make" ]; then
    cd "$HPX_BUILD_PATH"
    MAKECMD="make -k -j$NUM_BUILDCORES"
    INSTALLCMD="make install"
    
    $MAKECMD
    $INSTALLCMD
fi
