#!/bin/bash

if [ -z ${DGSWEMV2_ROOT+x} ]; then
    DGSWEMV2_ROOT_="${HOME}/dgswemv2"
else
    DGSWEMV2_ROOT_=$DGSWEMV2_ROOT
fi

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
make MANUFACTURED_SOLUTION_SERIAL
make MANUFACTURED_SOLUTION_HPX
make MANUFACTURED_SOLUTION_OMPI
echo ""
echo "Setting up runtime files..."
cd $HOME
mkdir -p dgswemv2_test
cp -r $DGSWEMV2_ROOT_/examples/manufactured_solution/input_files/* dgswemv2_test

cd dgswemv2_test
#Halve the manufactured solution run time to shorten circleci test time
sed -i 's/  end_time: 3600                #in seconds/  end_time: 1800/g' dgswemv2_input.15
$DGSWEMV2_ROOT_/build/mesh_generators/rectangular_mesh_generator mesh_generator_input.yml

echo "Running Serial Test case..."
rm -f serial.out
$DGSWEMV2_ROOT_/build/examples/MANUFACTURED_SOLUTION_SERIAL dgswemv2_input.15 &> serial.out

echo "Running HPX Test case..."
rm -f hpx.out
$DGSWEMV2_ROOT_/build/partitioner/partitioner dgswemv2_input.15 2 1
$DGSWEMV2_ROOT_/build/examples/MANUFACTURED_SOLUTION_HPX dgswemv2_input_parallelized.15 --hpx:threads=2 &> hpx.out

echo "Running MPI Test case..."
rm rectangular_mesh_*
rm -f ompi.out
$DGSWEMV2_ROOT_/build/partitioner/partitioner dgswemv2_input.15 2 1 2
#Since OpenMPI does not by default support MPI_THREAD_MULTIPLE
# we set OMP_NUM_THREADS=1
OMP_NUM_THREADS=1 mpirun -np 2 $DGSWEMV2_ROOT_/build/examples/MANUFACTURED_SOLUTION_OMPI dgswemv2_input_parallelized.15 &> ompi.out

python $DGSWEMV2_ROOT_/scripts/correctness/compare_l2_errors.py
exit $?