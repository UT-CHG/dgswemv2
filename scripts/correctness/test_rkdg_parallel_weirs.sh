#!/bin/bash

if [ -z ${DGSWEMV2_ROOT+x} ]; then
    DGSWEMV2_ROOT_="${HOME}/dgswemv2"
else
    DGSWEMV2_ROOT_=$DGSWEMV2_ROOT
fi
DGSWEMV2_TEST="${HOME}/dgswemv2_parallel_weirs_test"

if [ $# -gt 1 ]; then
    echo "rkdg_parallel_weirs only accepts one optional parameter"
    echo "  location of the build directory relative to $DGSWEMV2_ROOT"
    echo "  the default is build"
    return 1
fi

if [ $# -eq 1 ]; then
    ABS_BUILD_DIR=$DGSWEMV2_ROOT_/${1}
else
    ABS_BUILD_DIR=$DGSWEMV2_ROOT_/build
fi


#exit the script if any command returns with non-zero status
set -e

echo "Running test to check that distributed weirs are correctly parallelized"
echo "Compiling code (if necessary)..."
cd $ABS_BUILD_DIR
make partitioner
num_build_cores=$(( $(nproc) - 1))

echo ""
echo "Setting up runtime files..."
mkdir -p ${DGSWEMV2_TEST}
cp -r $DGSWEMV2_ROOT_/test/files_for_testing/weir/* ${DGSWEMV2_TEST}

cd ${DGSWEMV2_TEST}

echo ""
echo "seting up serial test case"
MAIN_DIR="${DGSWEMV2_ROOT_}/source/"
sed -i.tmp '/return 0/i\
        simulation->ComputeL2Residual();\
' ${MAIN_DIR}/dgswemv2-serial.cpp
cd $ABS_BUILD_DIR
make -j ${num_build_cores} dgswemv2-serial
cd $MAIN_DIR
#undo the addition of the ComputeL2Residual call
mv dgswemv2-serial.cpp.tmp dgswemv2-serial.cpp
echo ""
echo "Running Serial Test case..."
cd $DGSWEMV2_TEST
rm -f serial.out
$ABS_BUILD_DIR/source/dgswemv2-serial dgswemv2_input.15 &> serial.out

echo ""
echo "Building HPX Test case..."
cd $DGSWEMV2_TEST
rm -f weir_*
$ABS_BUILD_DIR/partitioner/partitioner dgswemv2_input.15 3 1
sed -i.tmp '/    return hpx::finalize();/i\
    hpx::future<double> globalResidualL2 = ComputeL2Residual(simulation_clients);\
    std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(globalResidualL2.get()) << std::endl;\
' ${MAIN_DIR}/dgswemv2-hpx.cpp
cd $ABS_BUILD_DIR
make -j ${num_build_cores} dgswemv2-hpx
cd $MAIN_DIR
mv dgswemv2-hpx.cpp.tmp dgswemv2-hpx.cpp
echo ""
echo "Running HPX Test case..."
cd $DGSWEMV2_TEST
rm -f hpx.out
$ABS_BUILD_DIR/source/dgswemv2-hpx dgswemv2_input_parallelized.15 --hpx:threads=3 &> hpx.out

echo ""
echo "Building MPI Test case..."
sed -i.tmp '/        MPI_Finalize();/i\
        simulation->ComputeL2Residual();\
' ${MAIN_DIR}/dgswemv2-ompi.cpp
cd $ABS_BUILD_DIR
make -j ${num_build_cores} dgswemv2-ompi
cd $MAIN_DIR
mv dgswemv2-ompi.cpp.tmp dgswemv2-ompi.cpp

echo ""
echo "Running OMPI Test case..."
cd $DGSWEMV2_TEST
rm -f weir_*
$ABS_BUILD_DIR/partitioner/partitioner dgswemv2_input.15 2 1 2
rm -f ompi.out
#Since OpenMPI does not by default support MPI_THREAD_MULTIPLE
# we set OMP_NUM_THREADS=1
#
# CI_MPI_CLI flag is an environemnt variable to run the mpi code as a root user.
# Running MPI as root is strongly discouraged by the OpemMPI people, so it really
# should only be used inside a container.
# See: https://github.com/open-mpi/ompi/issues/4451
OMP_NUM_THREADS=1 mpirun -np 2 ${CI_MPI_CLI} $ABS_BUILD_DIR/source/dgswemv2-ompi dgswemv2_input_parallelized.15 &> ompi.out

python $DGSWEMV2_ROOT_/scripts/correctness/compare_l2_errors.py
exit $?