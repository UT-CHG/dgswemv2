#include "utilities/linear_algebra/vector_dynamic.hpp"
#include "utilities/linear_algebra/matrix_dynamic.hpp"
#include "utilities/almost_equal.hpp"
#include "general_definitions.hpp"

int main() {
    bool error_found = false;

    using Utilities::almost_equal;

    Utilities::LinearAlgebra::VectorDyn<double> vec(2);
    Utilities::LinearAlgebra::VectorDyn<double> vec_res(3);

    vec(0) = 1;
    vec(1) = 2;

    // this should work if all overloads defined properly and result in -2*vec
    vec_res = vec - (vec + (vec + vec) + vec + 1.0 * vec + vec * 1.0 - vec - (vec - vec) - vec - vec / 1.0);

    if (!almost_equal(vec_res(0), -2.0, 100) || !almost_equal(vec_res(1), -4.0, 100)) {
        error_found = true;
    }

    // this should work if all overloads defined properly and result in -2*vec
    vec_res = -vec / 1.0;
    vec_res -= vec;
    vec_res -= vec - vec;
    vec_res -= vec;
    vec_res += vec * 1.0;
    vec_res += 1.0 * vec;
    vec_res += vec;
    vec_res += vec + vec;
    vec_res += vec;
    vec_res *= -1.0;
    vec_res += vec;
    vec_res /= 1.0;

    if (!almost_equal(vec_res(0), -2.0, 100) || !almost_equal(vec_res(1), -4.0, 100)) {
        error_found = true;
    }

    // check transpose, expr * vec , vec * expr
    Utilities::LinearAlgebra::MatrixDyn<double> outer_product(2, 2);
    outer_product = vec * vec_res.transpose();

    if (!almost_equal(outer_product(0, 0), -2.0, 100) || !almost_equal(outer_product(0, 1), -4.0, 100) ||
        !almost_equal(outer_product(1, 0), -4.0, 100) || !almost_equal(outer_product(1, 1), -8.0, 100)) {
        error_found = true;
    }

    double inner_product = vec.transpose() * vec_res;

    if (!almost_equal(inner_product, -10.0, 100)) {
        error_found = true;
    }

    vec_res = vec;  // check assignment

    if (!almost_equal(vec_res(0), 1.0, 100) || !almost_equal(vec_res(1), 2.0, 100)) {
        error_found = true;
    }

    double dot = vec.dot(vec_res);  // check dot

    if (!almost_equal(dot, 5.0, 100)) {
        error_found = true;
    }

    Utilities::LinearAlgebra::MatrixDyn<double> mat(2, 2);
    Utilities::LinearAlgebra::MatrixDyn<double> mat_res(2, 2);

    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(1, 0) = 3;
    mat(1, 1) = 4;

    // this should work if all overloads defined properly and result in -2*mat
    mat_res = mat - (mat + (mat + mat) + mat + 1.0 * mat + mat * 1.0 - mat - (mat - mat) - mat - mat / 1.0);

    if (!almost_equal(mat_res(0, 0), -2.0, 100) || !almost_equal(mat_res(0, 1), -4.0, 100) ||
        !almost_equal(mat_res(1, 0), -6.0, 100) || !almost_equal(mat_res(1, 1), -8.0, 100)) {
        error_found = true;
    }

    // this should work if all overloads defined properly and result in -2*vec
    mat_res = -mat / 1.0;
    mat_res -= mat;
    mat_res -= mat - mat;
    mat_res -= mat;
    mat_res += mat * 1.0;
    mat_res += 1.0 * mat;
    mat_res += mat;
    mat_res += mat + mat;
    mat_res += mat;
    mat_res *= -1.0;
    mat_res += mat;
    mat_res /= 1.0;

    if (!almost_equal(mat_res(0, 0), -2.0, 100) || !almost_equal(mat_res(0, 1), -4.0, 100) ||
        !almost_equal(mat_res(1, 0), -6.0, 100) || !almost_equal(mat_res(1, 1), -8.0, 100)) {
        error_found = true;
    }

    mat_res = mat;  // check assignment

    if (!almost_equal(mat_res(0, 0), 1.0, 100) || !almost_equal(mat_res(0, 1), 2.0, 100) ||
        !almost_equal(mat_res(1, 0), 3.0, 100) || !almost_equal(mat_res(1, 1), 4.0, 100)) {
        error_found = true;
    }

    mat_res = mat.transpose();  // check transpose

    if (!almost_equal(mat_res(0, 0), 1.0, 100) || !almost_equal(mat_res(0, 1), 3.0, 100) ||
        !almost_equal(mat_res(1, 0), 2.0, 100) || !almost_equal(mat_res(1, 1), 4.0, 100)) {
        error_found = true;
    }

    // check mat * mat
    mat_res = mat * mat;

    if (!almost_equal(mat_res(0, 0), 7.0, 100) || !almost_equal(mat_res(0, 1), 10.0, 100) ||
        !almost_equal(mat_res(1, 0), 15.0, 100) || !almost_equal(mat_res(1, 1), 22.0, 100)) {
        error_found = true;
    }

    // check mat * vec
    vec_res = mat * vec;

    if (!almost_equal(vec_res(0), 5.0, 100) || !almost_equal(vec_res(1), 11.0, 100)) {
        error_found = true;
    }

    // check vec^T * mat
    vec_res = vec.transpose() * mat;

    if (!almost_equal(vec_res(0), 7.0, 100) || !almost_equal(vec_res(1), 10.0, 100)) {
        error_found = true;
    }

    // check vec * mat (outer product)
    Utilities::LinearAlgebra::MatrixDyn<double> mat_out(1, 2);

    mat_out(0, 0) = -2;
    mat_out(0, 1) = -4;

    outer_product = vec * mat_out;

    if (!almost_equal(outer_product(0, 0), -2.0, 100) || !almost_equal(outer_product(0, 1), -4.0, 100) ||
        !almost_equal(outer_product(1, 0), -4.0, 100) || !almost_equal(outer_product(1, 1), -8.0, 100)) {
        error_found = true;
    }

    if (error_found) {
        return 1;
    }
    return 0;
};
