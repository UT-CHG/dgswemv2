#!/bin/bash

DGSWEMV2_ROOT="${HOME}/dgswemv2"
DGSWEMV2_TEST="${HOME}/dgswemv2_parallel_weirs_test"

echo "Running test to check that distributed weirs are correctly parallelized"
echo "Compiling code (if necessary)..."
cd $DGSWEMV2_ROOT/build
make partitioner
#make DG_HYPER_SWE_HPX
#make DG_HYPER_SWE_OMPI
echo ""
echo "Setting up runtime files..."
mkdir -p ${DGSWEMV2_TEST}
cp -r $DGSWEMV2_ROOT/test/files_for_testing/weir/* ${DGSWEMV2_TEST}

cd ${DGSWEMV2_TEST}
#Halve the manufactured solution run time to shorten circleci test time
#sed -i 's/  end_time: 3600                #in seconds/  end_time: 1800/g' dgswemv2_input.15

echo ""
echo "seting up serial test case"
MAIN_DIR="${DGSWEMV2_ROOT}/source/problem/SWE/main_files"
sed -i.tmp '/return 0/i\
        simulation.ComputeL2Residual();\
' ${MAIN_DIR}/main_swe.cpp
cd ${DGSWEMV2_ROOT}/build
make DG_HYPER_SWE_SERIAL
cd $MAIN_DIR
#undo the addition of the ComputeL2Residual call
mv main_swe.cpp.tmp main_swe.cpp
echo ""
echo "Running Serial Test case..."
cd $DGSWEMV2_TEST
rm -f serial.out
$DGSWEMV2_ROOT/build/source/DG_HYPER_SWE_SERIAL dgswemv2_input.15 &> serial.out

echo ""
echo "Building HPX Test case..."
cd $DGSWEMV2_TEST
rm weir_*
$DGSWEMV2_ROOT/build/partitioner/partitioner dgswemv2_input.15 3 1
sed -i.tmp '/    return hpx::finalize();/i\
    hpx::future<double> globalResidualL2 = ComputeL2Residual(simulation_clients);\
    std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(globalResidualL2.get()) << std::endl;\
' ${MAIN_DIR}/hpx_main_swe.cpp
cd ${DGSWEMV2_ROOT}/build
make DG_HYPER_SWE_HPX
cd $MAIN_DIR
mv hpx_main_swe.cpp.tmp hpx_main_swe.cpp
echo ""
echo "Running HPX Test case..."
cd $DGSWEMV2_TEST
rm -f hpx.out
$DGSWEMV2_ROOT/build/source/DG_HYPER_SWE_HPX dgswemv2_input_parallelized.15 --hpx:threads=3 &> hpx.out

echo ""
echo "Building MPI Test case..."
sed -i.tmp '/        MPI_Finalize();/i\
        simulation.ComputeL2Residual();\
' ${MAIN_DIR}/ompi_main_swe.cpp
cd ${DGSWEMV2_ROOT}/build
make DG_HYPER_SWE_OMPI
cd $MAIN_DIR
mv ompi_main_swe.cpp.tmp ompi_main_swe.cpp



echo ""
echo "Running OMPI Test case..."
cd $DGSWEMV2_TEST
rm weir_*
$DGSWEMV2_ROOT/build/partitioner/partitioner dgswemv2_input.15 3 1 3
rm -f ompi.out
#Since OpenMPI does not by default support MPI_THREAD_MULTIPLE
# we set OMP_NUM_THREADS=1
OMP_NUM_THREADS=1 mpirun -np 3 $DGSWEMV2_ROOT/build/source/DG_HYPER_SWE_OMPI dgswemv2_input_parallelized.15 &> ompi.out

python $DGSWEMV2_ROOT/scripts/correctness/compare_l2_errors.py
exit $?