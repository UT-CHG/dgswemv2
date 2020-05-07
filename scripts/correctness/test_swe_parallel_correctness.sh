#!/bin/bash

if [ -z ${DGSWEMV2_ROOT+x} ]; then
    DGSWEMV2_ROOT_="${HOME}/dgswemv2"
else
    DGSWEMV2_ROOT_=$DGSWEMV2_ROOT
fi

if [ $# -gt 2 ]; then
    echo "parallel_correctness_only accepts two optional parameters"
    echo "  (1) problem type"
    echo "      (default: rkdg_swe)"
    echo "  (2) location of the build directory relative to $DGSWEMV2_ROOT"
    echo "      (default: build)"

    return 1
fi

if [ $# -gt 0 ]; then
    PROBLEM=${1}
else
    PROBLEM="rkdg_swe"
fi

if [ $# -eq 2 ]; then
    BUILD_DIR=${2}
else
    BUILD_DIR="build"
fi


#exit the script if any command returns with non-zero status
set -e

source ${DGSWEMV2_ROOT_}/examples/swe_manufactured_solution/build.sh

echo "Building mesh for manufactured solution..."
cd $DGSWEMV2_ROOT_/mesh_generators
cat > bathymetry.hpp <<EOL
#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    return 2.0;  // bathymetry for manufactured solution
}
#endif
EOL
echo ""
echo "Compiling code (if necessary)..."
cd $DGSWEMV2_ROOT_/${BUILD_DIR}
make quad_mesh_generator
make partitioner
make_swe_manufactured_solution ${DGSWEMV2_ROOT_} serial ${BUILD_DIR}
make_swe_manufactured_solution ${DGSWEMV2_ROOT_} ompi ${BUILD_DIR}
echo ""
echo "Setting up runtime files..."
cd $HOME
mkdir -p dgswemv2_test
cp -r $DGSWEMV2_ROOT_/examples/swe_manufactured_solution/input_files/* dgswemv2_test

cd dgswemv2_test
sed -i "s/  name: rkdg_swe/  name: ${PROBLEM}/g" dgswemv2_input.15
$DGSWEMV2_ROOT_/$BUILD_DIR/mesh_generators/quad_mesh_generator mesh_generator_input.yml

echo "Running Serial Test case..."
rm -f serial.out
$DGSWEMV2_ROOT_/$BUILD_DIR/source/manufactured-solution-swe-serial dgswemv2_input.15 &> serial.out

echo "Running MPI Test case..."
rm -f ompi.out
$DGSWEMV2_ROOT_/$BUILD_DIR/partitioner/partitioner dgswemv2_input.15 2 1 2
#Since OpenMPI does not by default support MPI_THREAD_MULTIPLE
# we set OMP_NUM_THREADS=1
#
# CI_MPI_CLI flag is an environemnt variable to run the mpi code as a root user.
# Running MPI as root is strongly discouraged by the OpemMPI people, so it really
# should only be used inside a container.
# See: https://github.com/open-mpi/ompi/issues/4451
OMP_NUM_THREADS=1 mpirun -np 2 ${CI_MPI_CLI} $DGSWEMV2_ROOT_/$BUILD_DIR/source/manufactured-solution-swe-ompi dgswemv2_input_parallelized.15 &> ompi.out

python $DGSWEMV2_ROOT_/scripts/correctness/compare_l2_errors.py
exit $?