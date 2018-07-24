#ifndef USE_EIGEN_HPP
#define USE_EIGEN_HPP

#include "general_definitions.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

template <typename T, uint m>
using StatVector = Eigen::Matrix<T, m, 1>;
template <typename T, uint m, uint n>
using StatMatrix = Eigen::Matrix<T, m, n>;

template <typename T>
using DynVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using DynRowVector = Eigen::Matrix<T, 1, Eigen::Dynamic>;
template <typename T>
using DynMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T, uint m>
using HybMatrix = Eigen::Matrix<T, m, Eigen::Dynamic>;

template <typename T>
using SparseVector = Eigen::SparseVector<T>;
template <typename T>
using SparseMatrix = Eigen::SparseMatrix<T>;

template <typename T>
decltype(auto) IdentityMatrix(uint size) {
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(size, size);
}

template <typename T>
DynVector<T> IdentityVector(uint size) {
    DynVector<double> I_vector = DynVector<double>::Zero(size * size);

    for (uint i = 0; i < size; i++) {
        I_vector[i * size + i] = 1.0;
    }

    return I_vector;
}

/* Vector/Matrix (aka Tensor) Operations */
template <typename ArrayType>
void set_constant(ArrayType& array, double value) {
    array = ArrayType::Constant(array.rows(), array.cols(), value);
}

template <typename ArrayType>
void set_constant(ArrayType&& array, double value) {
    array = ArrayType::Constant(array.rows(), array.cols(), value);
}

template <typename ArrayType>
decltype(auto) transpose(const ArrayType& array) {
    return array.transpose();
}

template <typename ArrayType>
double norm(const ArrayType& array) {
    return array.norm();
}

template <typename LeftArrayType, typename RightArrayType>
decltype(auto) cwise_multiplication(const LeftArrayType& array_left, const RightArrayType& array_right) {
    return array_left.cwiseProduct(array_right);
}

template <typename LeftArrayType, typename RightArrayType>
decltype(auto) cwise_division(const LeftArrayType& array_left, const RightArrayType& array_right) {
    return array_left.cwiseQuotient(array_right);
}

/* Vector Operations */
template <typename VectorType>
decltype(auto) subvector(VectorType& vector, uint start_row, uint size_row) {
    return vector.segment(start_row, size_row);
}

template <typename T, uint n>
DynMatrix<T> reshape_jacobian_vector(const DynVector<T>& vector) {
    return DynMatrix<T>(n, n)();
}

/* Matrix Operations */
template <typename MatrixType>
uint rows(const MatrixType& matrix) {
    return matrix.rows();
}

template <typename MatrixType>
uint columns(const MatrixType& matrix) {
    return matrix.cols();
}

template <typename MatrixType>
decltype(auto) submatrix(MatrixType& matrix, uint start_row, uint start_col, uint size_row, uint size_col) {
    return matrix.block(start_row, start_col, size_row, size_col);
}

template <typename MatrixType>
decltype(auto) row(MatrixType& matrix, uint row) {
    return matrix.row(row);
}

template <typename MatrixType>
decltype(auto) column(MatrixType& matrix, uint col) {
    return matrix.col(col);
}

template <typename MatrixType>
decltype(auto) determinant(MatrixType& matrix) {
    return matrix.determinant();
}

template <typename MatrixType>
decltype(auto) inverse(MatrixType& matrix) {
    return matrix.inverse();
}

/* Solving Linear System */
template <typename MatrixType, typename VectorType>
void solve_sle(MatrixType& A, VectorType& b) {
    b = A.fullPivLu().solve(b);
}

template <typename VectorType, typename T>
void solve_sle(SparseMatrix<T>& A_sparse, VectorType& b) {}

#endif