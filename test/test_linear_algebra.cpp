#include "utilities/linear_algebra.hpp"

#include <iostream>
//#include <petscksp.h>

int main(int argc, char* args[]) {
    DynVector<double> a(9);
    DynVector<double> a2(9);
#ifdef USE_BLAZE
    a  = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    a2 = {5, 5, 5, 5, 5, 5, 5, 5, 5};

    DynMatrix<double> A = reshape<double, SO::ColumnMajor>(a, 3, 3);

    DynVector<double> b = flatten<double>(A);

    std::cout << A << a << b << "works\n";


    DynVector<double> c(4);
    c = { 1,2,3,4};
    std::cout << "         b: " << c << std::endl;
    std::cout << "reverse(b): " << reverse(c) << std::endl;

    DynVector<double> d(4);
    d = {2,1,2,1};
    std::cout << " d: " << d << std::endl;
    std::cout << " d + reverse(c): " << d + reverse(c) << std::endl;
    std::cout << " reverse(c + d): " << reverse(c+d) << std::endl;

    auto e = c + d;
    static_assert( std::is_same<decltype(e)::ResultType, DynVector<double>>::value);
#endif

#ifdef USE_EIGEN
    a  << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    a2 << 5, 5, 5, 5, 5, 5, 5, 5, 5;
    DynMatrix<double> A = reshape<double, 3>(a);

    std::cout << a << '\n' << A << "\nworks\n";
#endif

    bool error_found{false};
    { // check make rows functionality
        std::array<DynMatrix<double>,2> As{A,A};
        std::array<DynView<double>,2> Arows = make_rows_as_views(As,2);

        for ( uint i = 0; i < columns(A); ++i ) {
            if ( Arows[0][i] != Arows[1][i] ) {
                error_found = true;
                std::cout << "At row 2, index " << i << '\n'
                          << " Arows[0][i] = " << Arows[0][i] << '\n'
                          << " Arows[1][i] = " << Arows[1][i] << '\n';
            }
        }
    }

    { //check that row functionality works properly
        DynRow<double> r = row(A,2);
        A(2,2) = 600;
        if ( A(2,2) != r[2] ) {
            error_found = true;
            std::cout << "row not holding reference to A's data\n"
                      << " row(A,2)[2] = " << r[2] << "(should be " << A(2,2) << ")\n";
        }

    }

    { //check max_vec works properly
        DynVector<double> a_max = max_vec(a,a2);
        for ( uint i = 0; i < a_max.size(); ++i) {
            if ( a_max[i] < std::max(a[i], a2[i]) ) {
                error_found = true;
                std::cout << "Error found in max_vec with 2 arguments\n"
                          << "  " << a_max[i] << " < max( " << a[i] << ", " << a2[i] << ")\n";
            }
        }

        double a_maxval = max_vec(a);
        if ( a_maxval != 9 ) {
            error_found = true;
            std::cout << "Error: max value in a is 9, not " << a_maxval << std::endl;
        }
    }

    { //check for component-wise matrix division
        DynMatrix<double> A2 = reshape<double, SO::ColumnMajor>(a2, 3, 3);
        std::cout << " a/a2 = " << mat_cw_div(A,A2) << std::endl;
    }
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

    return error_found;
}
