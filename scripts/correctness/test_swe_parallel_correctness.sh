#!/bin/bash

if [ -z ${DGSWEMV2_ROOT+x} ]; then
    DGSWEMV2_ROOT_="${HOME}/dgswemv2"
else
    DGSWEMV2_ROOT_=$DGSWEMV2_ROOT
fi

if [ $# -gt 1 ]; then
    echo "parallel_correctness_only accepts one optional parameter"
    echo "  specifying the problem type"
    return 1
fi

if [ $# -eq 1 ]; then
    PROBLEM=${1}
else
    PROBLEM="rkdg_swe"
fi

source $DGSWEMV2_ROOT/examples/swe_manufactured_solution/build.sh

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
cd $DGSWEMV2_ROOT_/build
make rectangular_mesh_generator
make partitioner
make_swe_manufactured_solution ${DGSWEMV2_ROOT} serial
make_swe_manufactured_solution ${DGSWEMV2_ROOT} ompi
make_swe_manufactured_solution ${DGSWEMV2_ROOT} hpx
echo ""
echo "Setting up runtime files..."
cd $HOME
mkdir -p dgswemv2_test
cp -r $DGSWEMV2_ROOT_/examples/swe_manufactured_solution/input_files/* dgswemv2_test

cd dgswemv2_test
#Halve the manufactured solution run time to shorten circleci test time
sed -i 's/  end_time: 25-11-1987 01:00:00                /  end_time: 25-11-1987 00:30:00/g' dgswemv2_input.15
sed -i "s/  name: rkdg_swe/  name: ${PROBLEM}/g" dgswemv2_input.15
$DGSWEMV2_ROOT_/build/mesh_generators/rectangular_mesh_generator mesh_generator_input.yml

echo "Running Serial Test case..."
rm -f serial.out
$DGSWEMV2_ROOT_/build/source/manufactured-solution-swe-serial dgswemv2_input.15 &> serial.out

echo "Running HPX Test case..."
rm -f hpx.out
$DGSWEMV2_ROOT_/build/partitioner/partitioner dgswemv2_input.15 2 1
$DGSWEMV2_ROOT_/build/source/manufactured-solution-swe-hpx dgswemv2_input_parallelized.15 --hpx:threads=2 &> hpx.out

echo "Running MPI Test case..."
rm rectangular_mesh_*
rm -f ompi.out
$DGSWEMV2_ROOT_/build/partitioner/partitioner dgswemv2_input.15 2 1 2
#Since OpenMPI does not by default support MPI_THREAD_MULTIPLE
# we set OMP_NUM_THREADS=1
#
# CI_MPI_CLI flag is an environemnt variable to run the mpi code as a root user.
# Running MPI as root is strongly discouraged by the OpemMPI people, so it really
# should only be used inside a container.
# See: https://github.com/open-mpi/ompi/issues/4451
OMP_NUM_THREADS=1 mpirun -np 2 ${CI_MPI_CLI} $DGSWEMV2_ROOT_/build/source/manufactured-solution-swe-ompi dgswemv2_input_parallelized.15 &> ompi.out

python $DGSWEMV2_ROOT_/scripts/correctness/compare_l2_errors.py
exit $?