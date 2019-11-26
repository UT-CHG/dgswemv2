 #!/bin/bash

################################################################################
# print usage information
# Inputs:
#   1. script name
function usage {
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

################################################################################
# parses inputs to build scripts (looking for correctness)
# Inputs:
#  1. Simply pass all command line arguments with "$@"
# Outputs:
#  1. Sets CONFIGFILE
function parse_args {
    CONFIGFILE=""

    if [ "$#" -gt 2 ]; then usage; fi

    if [ "$#" == 2 ] && [ "$1" == "-c" ]; then
        CONFIGFILE=$2
        return
    fi

    if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then usage; fi

    # Fix me: handle more options
    if [ "$#" -gt 0 ]; then
        if [ "$1" != "clean" ] && [ "$1" != "no-make" ] && [ "$1" != "reinstall" ]; then
	    echo "invalid option: $1"
	    usage
        fi
    fi
}
################################################################################
# see if system supports modules and if so load modules specified in $MODULES
#Inputs:
# 1. Space delimited string of modules to be loaded
function try_loading_modules {
if [ -x "$(command -v modules list)" ]; then
    echo "MODULES = ${1}"

    module purge

    for module in ${1}; do
        module load $module
    done

    module list
fi
}

################################################################################
# load config file
#Inputs:
# 1. Config file
function load_config_file {

    echo -n "config file: "
    $START_BOLD
    echo "${1}"
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
    done < ${1}

    echo "Press enter to continue with these settings (or ctrl-c to exit)."
    read answer

    if [ -s "${1}" ]; then
	source ${1}
    else
	echo "config file ${1} either empty or non-existent"
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
}

################################################################################
# Clean up build directory
# Inputs:
#  1. Build directory
#  2. clean up command
function clean_up {
    $START_BOLD
    echo "  clean:"
    $END_BOLD
    echo "about to execute:"
    echo "rm -rf ${1}"
    $START_BOLD
    echo "is this okay? [y/N]"
    $END_BOLD
    read answer
    if echo "$answer" | grep -iq "^y" ;then
	echo "removing build directory ${1}"
	rm -rf "${1}"
    else
	echo "doing nothing."
    fi
    exit 0
}