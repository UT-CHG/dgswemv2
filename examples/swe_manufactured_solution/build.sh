#!/bin/bash

########################################
# Build script for manufactured solution
# Takes arguments:
#  $1: path/to/dgswemv2
#  $2: (optional) target type
#      must be one of (serial(default), ompi, hpx)
#  $3: (optional) build_directory
#      defaults to $1/build
make_swe_manufactured_solution() {
    if [ $# -lt 4 -a $# -gt 0 ]; then
	DGSWEMV2_DIR=${1}
	if [ $# -gt 1 ]; then
	    TARGET=${2}
	else
	    TARGET="serial"
	fi
	if [ $# -gt 2 ]; then
	    #TODO add logic to check for absolute paths
	    BUILD_DIR=${DGSWEMV2_DIR}/${3}
	else
	    BUILD_DIR=${DGSWEMV2_DIR}/build
	fi

	#check that directories exist
	if [ ! -d "$DGSWEMV2_DIR" ]; then
	    echo "Error: could not find dgswemv2 directory: ${DGSWEMV2_DIR}"
	    return 4
	fi

	if [ ! -d "$BUILD_DIR" ]; then
	    echo "Error: could not find build directory: ${BUILD_DIR}"
	    return 3
	fi

	#check that target is valid
	if [ "${TARGET}" != "serial" ] && \
	   [ "${TARGET}" != "ompi" ] && \
	   [ "${TARGET}" != "hpx" ] ; then
	    echo "Error: invalid target type: ${TARGET}"
	    echo "       please select one of serial ompi or hpx"
	    return 2
	fi

	#Print out input state
	echo "DGSWEMV2_DIR: ${DGSWEMV2_DIR}"
	echo "TARGET: ${TARGET}"
	echo "BUILD_DIR: ${BUILD_DIR}"

	#swap files
	MANSOL_DIR=${DGSWEMV2_DIR}/examples/swe_manufactured_solution
	SWE_DIR=${DGSWEMV2_DIR}/source/problem/SWE/problem_function_files

	cp ${SWE_DIR}/swe_initial_condition_functions.cpp ${SWE_DIR}/swe_initial_condition_functions.cpp.tmp
	cp ${SWE_DIR}/swe_true_solution_functions.cpp ${SWE_DIR}/swe_true_solution_functions.cpp.tmp
	cp ${SWE_DIR}/swe_source_functions.cpp ${SWE_DIR}/swe_source_functions.cpp.tmp

	cp ${MANSOL_DIR}/manufactured_swe_initial_condition_functions.cpp ${SWE_DIR}/swe_initial_condition_functions.cpp
	cp ${MANSOL_DIR}/manufactured_swe_true_solution_functions.cpp ${SWE_DIR}/swe_true_solution_functions.cpp
	cp ${MANSOL_DIR}/manufactured_swe_source_functions.cpp ${SWE_DIR}/swe_source_functions.cpp

	MAIN_DIR="${DGSWEMV2_DIR}/source/"

	if [ "${TARGET}" == "serial" ]; then
	    sed -i.tmp '/return 0/i\
	    simulation->ComputeL2Residual();' ${MAIN_DIR}/dgswemv2-serial.cpp
	elif [ "${TARGET}" == "ompi" ]; then
	    sed -i.tmp '/        MPI_Finalize();/i\
                simulation->ComputeL2Residual();\' ${MAIN_DIR}/dgswemv2-ompi.cpp
	elif [ "${TARGET}" == "hpx" ]; then
	    sed -i.tmp '/    return hpx::finalize();/i\
                hpx::future<double> globalResidualL2 = ComputeL2Residual(simulation_clients);\
                std::cout << "L2 error: " << std::setprecision(15) << std::sqrt(globalResidualL2.get()) << std::endl;\
                ' ${MAIN_DIR}/dgswemv2-hpx.cpp
	fi

	old_dir=${PWD}
	cd $BUILD_DIR
	make dgswemv2-${TARGET}
	status=$?
	mv source/dgswemv2-${TARGET} source/manufactured-solution-swe-${TARGET}
	cd $old_dir

	#cleaning up
	mv ${SWE_DIR}/swe_initial_condition_functions.cpp.tmp ${SWE_DIR}/swe_initial_condition_functions.cpp
	mv ${SWE_DIR}/swe_true_solution_functions.cpp.tmp ${SWE_DIR}/swe_true_solution_functions.cpp
	mv ${SWE_DIR}/swe_source_functions.cpp.tmp ${SWE_DIR}/swe_source_functions.cpp
	mv ${MAIN_DIR}/dgswemv2-${TARGET}.cpp.tmp ${MAIN_DIR}/dgswemv2-${TARGET}.cpp

	return ${status}
    else
	echo 'Usage:'
	echo '    Build script for manufactured solution'
	echo '    Takes arguments:'
	echo '      ${1}: path/to/dgswemv2'
	echo '      ${2}: (optional) target type'
	echo '            must be one of (serial(default), ompi, hpx)'
	echo '      ${3}: (optional) build_directory'
	echo '            defaults to ${1}/build'
	return 1
    fi

}