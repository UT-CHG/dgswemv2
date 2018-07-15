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
using DynMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using SparseVector = Eigen::SparseVector<T>;
template <typename T>
using SparseMatrix = Eigen::SparseMatrix<T>;

template <typename T>
decltype(auto) IdentityMatrix(uint size) {
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(size, size);
}

/* Vector/Matrix (aka Tensor) Operations */
template <typename TensorType>
double norm(TensorType& tensor) {
    return tensor.norm();
}

template <typename TensorType>
double sqr_norm(TensorType& tensor) {
    return tensor.squaredNorm();
}

/* Vector Operations */
template <typename VectorType>
decltype(auto) subvector(VectorType& vector, uint start_row, uint size_row) {
    return vector.segment(start_row, size_row);
}

/* Matrix Operations */
template <typename MatrixType>
decltype(auto) submatrix(MatrixType& matrix, uint start_row, uint start_col, uint size_row, uint size_col) {
    return matrix.block(start_row, start_col, size_row, size_col);
}

template <typename MatrixType, uint col>
decltype(auto) column(MatrixType& matrix) {
    return matrix.col(col);
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