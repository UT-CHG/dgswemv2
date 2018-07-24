#ifndef USE_BLAZE_HPP
#define USE_BLAZE_HPP

#include "general_definitions.hpp"

#include <blaze/Math.h>
#include <blaze/math/Subvector.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/Column.h>
#include <blaze/math/Row.h>

template <typename T, uint m>
using StatVector = blaze::StaticVector<T, m>;
template <typename T, uint m, uint n>
using StatMatrix = blaze::StaticMatrix<T, m, n>;

template <typename T>
using DynVector = blaze::DynamicVector<T>;
template <typename T>
using DynRowVector = blaze::DynamicVector<T, blaze::rowVector>;
template <typename T>
using DynMatrix = blaze::DynamicMatrix<T>;
template <typename T, uint m>
using HybMatrix = blaze::HybridMatrix<T, m, 16>;

template <typename T>
using SparseVector = blaze::CompressedVector<T>;
template <typename T>
using SparseMatrix = blaze::CompressedMatrix<T>;

template <typename T>
using IdentityMatrix = blaze::IdentityMatrix<T>;

template <typename T>
DynVector<T> IdentityVector(uint size) {
    DynVector<double> I_vector(size * size, 0.0);

    for (uint i = 0; i < size; i++) {
        I_vector[i * size + i] = 1.0;
    }

    return I_vector;
}

/* Vector/Matrix (aka Array) Operations */
template <typename ArrayType>
void set_constant(ArrayType& array, double value) {
    array = value;
}

template <typename ArrayType>
void set_constant(ArrayType&& array, double value) {
    array = value;
}

template <typename ArrayType>
decltype(auto) transpose(const ArrayType& array) {
    return blaze::trans(array);
}

template <typename ArrayType>
double norm(const ArrayType& array) {
    return blaze::norm(array);
}

template <typename LeftArrayType, typename RightArrayType>
decltype(auto) cwise_multiplication(const LeftArrayType& array_left, const RightArrayType& array_right) {
    return array_left * array_right;
}

template <typename LeftArrayType, typename RightArrayType>
decltype(auto) cwise_division(const LeftArrayType& array_left, const RightArrayType& array_right) {
    return array_left / array_right;
}

/* Vector Operations */
template <typename VectorType>
decltype(auto) subvector(VectorType& vector, uint start_row, uint size_row) {
    return blaze::subvector(vector, start_row, size_row);
}

template <typename T, uint n>
DynMatrix<T> reshape_jacobian_vector(const DynVector<T>& vector) {
    return DynMatrix<T>(n, n, vector.data());
}

/* Matrix Operations */
template <typename MatrixType>
uint rows(const MatrixType& matrix) {
    return blaze::rows(matrix);
}

template <typename MatrixType>
uint columns(const MatrixType& matrix) {
    return blaze::columns(matrix);
}

template <typename MatrixType>
decltype(auto) submatrix(MatrixType& matrix, uint start_row, uint start_col, uint size_row, uint size_col) {
    return blaze::submatrix(matrix, start_row, start_col, size_row, size_col);
}

template <typename MatrixType>
decltype(auto) row(MatrixType& matrix, uint row) {
    return blaze::row(matrix, row);
}

template <typename MatrixType>
decltype(auto) column(MatrixType& matrix, uint col) {
    return blaze::column(matrix, col);
}

template <typename MatrixType>
double determinant(MatrixType& matrix) {
    return blaze::det(matrix);
}

template <typename MatrixType>
decltype(auto) inverse(MatrixType& matrix) {
    return blaze::inv(matrix);
}

/* Solving Linear System */
template <typename MatrixType, typename VectorType>
void solve_sle(MatrixType& A, VectorType& b) {
    int ipiv[b.size()];

    blaze::gesv(A, b, ipiv);
}

template <typename VectorType, typename T>
void solve_sle(SparseMatrix<T>& A_sparse, VectorType& b) {
    // Avoid using this function, use a library with sparse solvers, e.g. Eigen
    uint rows = A_sparse.rows();
    uint cols = A_sparse.columns();

    DynMatrix<double> A_dense(rows, cols, 0.0);

    A_dense = A_sparse;

    int ipiv[b.size()];

    blaze::gesv(A_dense, b, ipiv);
}

#endif