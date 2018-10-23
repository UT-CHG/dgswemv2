#include "general_definitions.hpp"
#include <petscksp.h>

/*

THIS TEST DOES NOT TEST ANYTHING
ITS HERE JUST TO CHECK FUNCTIONALITIES OF LINEAR ALGEBRA PACKAGES

*/

int main(int argc, char* args[]) {
#ifdef USE_BLAZE
    DynMatrix<double> A(3, 3);
    DynMatrix<double> B(3, 3);

    A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    B = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    DynMatrix<double> C = vec_cw_div(A, B);

    std::cout << A << B << C << "works\n";
#endif

#ifdef USE_EIGEN
    DynVector<double> a(9);
    a << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    DynMatrix<double> A = reshape<double, 3>(a);

    std::cout << a << '\n' << A << "\nworks\n";
#endif

    /*PetscInitialize(&argc, &args, (char*)0, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Vec vector;
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, 4, &vector);

    double values[4] = {1, 2, 3, 4};
    int indexes[4]   = {0, 1, 2, 3};

    if (rank == 0)
        VecSetValues(vector, 4, indexes, values, ADD_VALUES);

    VecAssemblyBegin(vector);
    VecAssemblyEnd(vector);

    VecView(vector, PETSC_VIEWER_STDOUT_WORLD);

    Vec vector_loc;
    VecCreateSeq(MPI_COMM_SELF, 6, &vector_loc);

    VecScatter scatter;
    IS from, to;

    int loc_indexes[6] = {0, 1, 2, 0, 1, 2};

    ISCreateGeneral(MPI_COMM_SELF, 6, loc_indexes, PETSC_COPY_VALUES, &from);
    ISCreateStride(MPI_COMM_SELF, 6, 0, 1, &to);

    VecScatterCreate(vector, from, vector_loc, to, &scatter);

    VecScatterBegin(scatter, vector, vector_loc, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, vector, vector_loc, INSERT_VALUES, SCATTER_FORWARD);

    ISDestroy(&from);
    ISDestroy(&to);
    VecScatterDestroy(&scatter);

    if (rank == 0)
        VecView(vector_loc, PETSC_VIEWER_STDOUT_SELF);

    Mat A;

    MatCreate(MPI_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 4, 4);

    MatSetUp(A);

    int indx_i[4] = {0, 1, 2, 3};
    int indx_j[4] = {0, 1, 2, 3};

    double vals[16] = {5, 6, 6, 8, 2, 2, 2, 8, 6, 6, 2, 8, 2, 3, 6, 7};

    if (rank == 0)
        MatSetValues(A, 4, indx_i, 4, indx_j, vals, INSERT_VALUES);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    Vec sol;
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, 4, &sol);

    KSP ksp;
    PC pc;

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);

    KSPSolve(ksp, vector, sol);

    VecView(sol, PETSC_VIEWER_STDOUT_WORLD);

    PetscFinalize();*/
}
