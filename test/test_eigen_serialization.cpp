#include <hpx/runtime/serialization/serialize.hpp>
#include "utilities/linear_algebra/use_eigen.hpp"
#include "utilities/almost_equal.hpp"

#include <iostream>

template <typename Matrix>
bool test_serialize_eigen(Matrix& i_mat) {
    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << i_mat;

    hpx::serialization::input_archive i_archive(buffer);
    Matrix o_mat;
    i_archive >> o_mat;

    bool error_found{false};

    Matrix residual = i_mat - o_mat;

    std::cout << "In object:\n" << i_mat << '\n'
              << "Out object:\n" << o_mat << '\n';

    error_found |= !Utilities::almost_equal(0., residual.norm());

    return error_found;
}

int main() {
    bool error_found{false};

    {  // serialize vector
        StatVector<double, 3> a;
        a << 1, 2, 3;
        error_found |= test_serialize_eigen(a);

        DynVector<double> b(3);
        b << 4, 5, 6;
        error_found |= test_serialize_eigen(b);
    }

    {  // serialize matrix

        StatMatrix<double, 3, 4> a;
        a << 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 7;
        error_found |= test_serialize_eigen(a);

        DynMatrix<double> b = DynMatrix<double>::Zero(4,4);
        b(0,0) = 1; b(0,1) = 2; b(0,2) = 3; b(0,3) = 4;
        b(2,0) = 1; b(2,1) = 2;
        b(3,0) = 1; b(3,1) = 2; b(3,2) = 3; b(3,3) = 5;
        error_found |= test_serialize_eigen(b);

        SparseMatrix<double> d(4,4);
        std::vector<Eigen::Triplet<double>> coefficients{{0,0,1},{0,1,2},{0,2,3},{0,3,4},{2,0,1},{2,2,2},{3,0,1},{3,1,2},{3,2,3},{3,3,5}};
        d.setFromTriplets(coefficients.begin(), coefficients.end());
        error_found |= test_serialize_eigen(d);
    }
    return error_found;
}