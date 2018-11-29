#ifndef USE_EIGEN_HPP
#define USE_EIGEN_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#define EIGEN_DONT_PARALLELIZE

#ifdef HAS_HPX
#include "serialization/eigen_matrix.hpp"
#endif

namespace SO {
constexpr int ColumnMajor = Eigen::StorageOptions::ColMajor;
constexpr int RowMajor    = Eigen::StorageOptions::RowMajor;
}


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

template <typename Matrix>
using Column = typename Eigen::DenseBase<Matrix>::ColXpr;

template <typename Matrix>
using Row = typename Eigen::DenseBase<Matrix>::RowXpr;

template <typename T>
using DynRow = Row<DynMatrix<T>>;

template <typename T>
using SparseVector = Eigen::SparseVector<T>;
template <typename T>
using SparseMatrix = Eigen::SparseMatrix<T>;

template <typename EigenType>
using AlignedAllocator = Eigen::aligned_allocator<EigenType>;

template <typename T>
DynMatrix<T> IdentityMatrix(const uint size) {
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(size, size);
}

template <typename T>
DynVector<T> IdentityVector(const uint size) {
    DynVector<T> I_vector = DynVector<T>::Zero(size * size);
    for (uint i = 0; i < size; ++i) {
        I_vector[i * size + i] = 1.0;
    }
    return I_vector;
}

template <typename T>
struct SparseMatrixMeta {
    std::vector<Eigen::Triplet<T>> data;

    void add_triplet(const uint row, const uint col, const T value) {
        this->data.emplace_back(Eigen::Triplet<T>(row, col, value));
    }

    void get_sparse_matrix(SparseMatrix<T>& sparse_matrix) { sparse_matrix.setFromTriplets(data.begin(), data.end()); }
};

/* Vector/Matrix (aka Tensor) Operations */
template <typename ArrayType>
void set_constant(ArrayType&& array, const double value) {
    array = std::remove_reference<ArrayType>::type::Constant(
        std::forward<ArrayType>(array).rows(), std::forward<ArrayType>(array).cols(), value);
}

template <typename ArrayType>
decltype(auto) transpose(const ArrayType& array) {
    return array.transpose();
}

template <typename ArrayType>
double norm(const ArrayType& array) {
    return array.norm();
}

template <typename ArrayType>
double sq_norm(const ArrayType& array) {
    return array.squaredNorm();
}

template <typename ArrayType>
decltype(auto) power(const ArrayType& array, const double exp) {
    return array.array().pow(exp);
}

/* Vector Operations */
template <typename LeftVectorType, typename RightVectorType>
decltype(auto) vec_cw_mult(const LeftVectorType& vector_left, const RightVectorType& vector_right) {
    return vector_left.cwiseProduct(vector_right);
}

template <typename LeftVectorType, typename RightVectorType>
decltype(auto) vec_cw_div(const LeftVectorType& vector_left, const RightVectorType& vector_right) {
    return vector_left.cwiseQuotient(vector_right);
}

template <typename T>
Eigen::Map<DynVector<T>> vector_from_array(T* array, const uint m) {
    return Eigen::Map<DynVector<T>>(array, m);
}

template <typename VectorType>
decltype(auto) subvector(VectorType&& vector, const uint start_row, const uint size_row) {
    return vector.segment(start_row, size_row);
}

template <typename T, int m, int n = m, int SO = Eigen::StorageOptions::RowMajor>
Eigen::Map<Eigen::Matrix<T, m, n, SO>> reshape(const StatVector<T, m * n>& vector) {
    return Eigen::Map<Eigen::Matrix<T, m, n, SO>>(const_cast<T*>(vector.data()), m, n);
}

template <typename T, int m, int SO = Eigen::StorageOptions::RowMajor>
Eigen::Map<Eigen::Matrix<T, m, Eigen::Dynamic, SO>> reshape(const DynVector<T>& vector, const int n) {
    return Eigen::Map<Eigen::Matrix<T, m, Eigen::Dynamic, SO>>(const_cast<T*>(vector.data()), m, n);
}

template <typename T, int SO = Eigen::StorageOptions::RowMajor>
Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, SO>> reshape(const DynVector<T>& vector,
                                                                         const int m,
                                                                         const int n) {
    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, SO>>(const_cast<T*>(vector.data()), m, n);
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
decltype(auto) submatrix(MatrixType&& matrix,
                         const uint start_row,
                         const uint start_col,
                         const uint size_row,
                         const uint size_col) {
    return matrix.block(start_row, start_col, size_row, size_col);
}

template <typename MatrixType>
decltype(auto) row(MatrixType&& matrix, const uint row) {
    return matrix.row(row);
}

template <typename MatrixType>
decltype(auto) column(MatrixType&& matrix, const uint col) {
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

template <typename T, int m, int n = m, int SO = Eigen::StorageOptions::RowMajor>
StatVector<T, m * n> flatten(const StatMatrix<T, m, n>& matrix) {
    StatVector<T, m * n> ret;
    Eigen::Map<Eigen::Matrix<double, m, n, SO>>(ret.data(), m, n) = matrix;
    return ret;
}

template <typename T, int m, int SO = Eigen::StorageOptions::RowMajor>
DynVector<T> flatten(const HybMatrix<T, m>& matrix) {
    uint n = matrix.cols();

    DynVector<T> ret(m * n);
    Eigen::Map<Eigen::Matrix<double, m, Eigen::Dynamic, SO>>(ret.data(), m, n) = matrix;
    return ret;
}

template <typename T, int SO = Eigen::StorageOptions::RowMajor>
DynVector<T> flatten(const DynMatrix<T>& matrix) {
    uint m = matrix.rows();
    uint n = matrix.cols();

    DynVector<T> ret(m * n);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, SO>>(ret.data(), m, n) = matrix;
    return ret;
}

/* Solving Linear System */
template <typename MatrixType, typename ArrayType>
void solve_sle(MatrixType& A, ArrayType& B) {
    B = A.fullPivLu().solve(B);
}

template <typename ArrayType, typename T>
void solve_sle(SparseMatrix<T>& A_sparse, ArrayType& B) {
    Eigen::SparseLU<SparseMatrix<T>> solver;
    solver.analyzePattern(A_sparse);
    solver.factorize(A_sparse);
    B = solver.solve(B);
}

#endif