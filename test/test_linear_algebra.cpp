#include "general_definitions.hpp"
#include <petscksp.h>

/*

THIS TEST DOES NOT TEST ANYTHING
ITS HERE JUST TO CHECK FUNCTIONALITIES OF LINEAR ALGEBRA PACKAGES

*/

int main(int argc, char* args[]) {
    /*#ifdef USE_BLAZE
        DynMatrix<double> A(3, 3);
        DynMatrix<double> B(3, 3);

        A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        B = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

        DynMatrix<double> C = cwise_division(A, B);

        std::cout << A << B << C << "works\n";
    #endif

    #ifdef USE_EIGEN
        DynVector<double> a(9);
        a << 1, 2, 3, 4, 5, 6, 7, 8, 9;

        DynMatrix<double> A = reshape<double, 3>(a);

        std::cout << a << '\n' << A << "\nworks\n";
    #endif*/

    PetscInitialize(&argc, &args, (char*)0, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Vec vector;
    VecCreateMPI(MPI_COMM_WORLD, 3, 12, &vector);

    double values[4] = {1.0, 2.0, 3.0, 4.0};
    int indexes[4]   = {2 * rank, 2 * rank + 1, 2 * rank + 2, 2 * rank + 3};

    // VecSet(vector, 1.0);
    VecSetValues(vector, 4, indexes, values, ADD_VALUES);

    VecAssemblyBegin(vector);
    VecAssemblyEnd(vector);

    VecView(vector, PETSC_VIEWER_STDOUT_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    PetscFinalize();
}
