#ifndef USE_BLAZE_HPP
#define USE_BLAZE_HPP

#include <blaze/Math.h>
#include <blaze/math/Subvector.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/Column.h>
#include <blaze/math/Row.h>
#include <blaze/util/AlignedAllocator.h>

#include <map>

#ifdef HAS_HPX
#include "serialization/blaze_vector.hpp"
#include "serialization/blaze_matrix.hpp"
#endif

using so_t = bool;

namespace SO {
constexpr bool ColumnMajor = blaze::columnMajor;
constexpr bool RowMajor    = blaze::rowMajor;
}

constexpr int hyb_mat_buff_size = 128;

template <typename T, uint m>
using StatVector = blaze::StaticVector<T, m>;
template <typename T, uint m, uint n>
using StatMatrix = blaze::StaticMatrix<T, m, n>;

template <typename T>
using DynVector = blaze::DynamicVector<T>;
template <typename T>
using DynRowVector = blaze::DynamicVector<T, blaze::rowVector>;
template <typename T, bool SO = blaze::rowMajor>
using DynMatrix = blaze::DynamicMatrix<T,SO>;
template <typename T, uint m>
using HybMatrix = blaze::HybridMatrix<T, m, hyb_mat_buff_size>;

template <typename Matrix>
using Column = decltype(blaze::column(std::declval<Matrix>(), std::declval<int>()));

template <typename Matrix>
using Row = decltype(blaze::row(std::declval<Matrix>(), std::declval<int>()));

template <typename T>
using DynRow = Row<DynMatrix<T>>;

template <typename Matrix>
using View = Row<Matrix>;

template <typename T, bool SO = blaze::rowMajor>
using DynView = View<DynMatrix<T, SO>>;

template <typename T>
using SparseVector = blaze::CompressedVector<T>;
template <typename T>
using SparseMatrix = blaze::CompressedMatrix<T>;

template <typename T>
using IdentityMatrix = blaze::IdentityMatrix<T>;

template <typename T, bool SO = blaze::rowMajor>
using DiagonalMatrix = blaze::DiagonalMatrix<blaze::CompressedMatrix<T, SO>>;

template <typename BlazeType>
using AlignedAllocator = blaze::AlignedAllocator<BlazeType>;

template <typename ArrayType>
struct Result {
    using type = typename ArrayType::ResultType;
};

template <typename T>
DynVector<T> IdentityVector(const uint size) {
    DynVector<T> I_vector(size * size, 0.0);
    for (uint i = 0; i < size; ++i) {
        I_vector[i * size + i] = 1.0;
    }
    return I_vector;
}

template <typename T>
struct SparseMatrixMeta {
    std::map<uint, std::map<uint, T>> data;

    void add_triplet(const uint row, const uint col, const T value) { this->data[row][col] = value; }

    void get_sparse_matrix(SparseMatrix<T>& sparse_matrix) {
        uint nel = 0;

        for (auto& row : data) {
            nel += row.second.size();
        }

        sparse_matrix.reserve(nel);

        for (uint row = 0; row < blaze::rows(sparse_matrix); ++row) {
            if (data.find(row) != data.end()) {
                for (auto& col : data.at(row)) {
                    sparse_matrix.append(row, col.first, col.second);
                }
                sparse_matrix.finalize(row);
            } else {
                sparse_matrix.finalize(row);
            }
        }
    }
};

/* Vector/Matrix (aka Array) Operations */
template <typename ArrayType>
void set_constant(ArrayType&& array, const double value) {
    array = value;
}

template <typename ArrayType>
decltype(auto) transpose(const ArrayType& array) {
    return blaze::trans(array);
}

template <typename ArrayType>
double sq_norm(const ArrayType& array) {
    return blaze::sqrNorm(array);
}

template <typename ArrayType>
double norm(const ArrayType& array) {
    return blaze::norm(array);
}

template <typename ArrayType, typename Number>
decltype(auto) pow_vec(const ArrayType& array, const Number exp) {
    return blaze::pow(array, exp);
}

template <typename ArrayType>
decltype(auto) sqrt_vec(const ArrayType& array) {
    return blaze::sqrt(array);
}

template <typename ArrayType>
decltype(auto) abs_vec(const ArrayType& array) {
    return blaze::abs(array);
}

template <typename ArrayType>
decltype(auto) cos_vec(const ArrayType& array) {
    return blaze::cos(array);
}

template <typename... ArrayTypes>
decltype(auto) max_vec(const ArrayTypes& ...arrays) {
    return blaze::max(arrays...);
}

/* Vector Operations */
template <typename LeftVectorType, typename RightVectorType>
decltype(auto) vec_cw_mult(const LeftVectorType& vector_left, const RightVectorType& vector_right) {
    return vector_left * vector_right;
}

template <typename LeftVectorType, typename RightVectorType>
decltype(auto) vec_cw_div(const LeftVectorType& vector_left, const RightVectorType& vector_right) {
    return vector_left / vector_right;
}

template <typename T>
DynVector<T> vector_from_array(T* array, const uint m) {
    return DynVector<T>(m, array);
}

template <typename VectorType>
decltype(auto) subvector(VectorType&& vector, const uint start_row, const uint size_row) {
    return blaze::subvector(std::forward<VectorType>(vector), start_row, size_row);
}

template <typename T, int m, int n = m, bool SO = blaze::rowMajor>
blaze::StaticMatrix<T, m, n, SO> reshape(const StatVector<T, m * n>& vector) {
    return blaze::StaticMatrix<T, m, n, SO>(m, n, vector.data());
}

template <typename T, int m, bool SO = blaze::rowMajor>
blaze::HybridMatrix<T, m, hyb_mat_buff_size, SO> reshape(const DynVector<T>& vector, const int n) {
    return blaze::HybridMatrix<T, m, hyb_mat_buff_size, SO>(m, n, vector.data());
}

template <typename T, bool SO = blaze::rowMajor>
blaze::DynamicMatrix<T, SO> reshape(const DynVector<T>& vector, const int m, const int n) {
    return blaze::DynamicMatrix<T, SO>(m, n, vector.data());
}

template <typename InputArrayType>
decltype(auto) reverse(InputArrayType& vector) {
    std::vector<size_t> reversed_indices(vector.size());
    for ( uint i= 0; i < vector.size(); ++i ) {
        reversed_indices[i] = vector.size()-1 - i;
    }
    return blaze::elements(vector, reversed_indices);
}

/* Matrix Operations */
template <typename T, bool SO = blaze::rowMajor>
blaze::DiagonalMatrix<blaze::CompressedMatrix<T, SO>> initialize_as_constant(size_t rows, T value) {
    //note that a diagonal matrix is always square
    blaze::DiagonalMatrix<blaze::CompressedMatrix<T, SO>> diag(rows, rows);
    for ( size_t i = 0UL; i < rows; ++i ) {
        diag.append(i, i, value);
        diag.finalize(i);
    }
    return diag;
}

template <typename MatrixType>
uint rows(const MatrixType& matrix) {
    return blaze::rows(matrix);
}

template <typename MatrixType>
uint columns(const MatrixType& matrix) {
    return blaze::columns(matrix);
}

template <typename MatrixType>
decltype(auto) submatrix(MatrixType&& matrix,
                         const uint start_row,
                         const uint start_col,
                         const uint size_row,
                         const uint size_col) {
    return blaze::submatrix(std::forward<MatrixType>(matrix), start_row, start_col, size_row, size_col);
}

template <typename MatrixType>
decltype(auto) row(MatrixType&& matrix, const uint row) {
    return blaze::row(std::forward<MatrixType>(matrix), row);
}

template <typename MatrixType>
decltype(auto) column(MatrixType&& matrix, const uint col) {
    return blaze::column(std::forward<MatrixType>(matrix), col);
}

template <typename MatrixType>
decltype(auto) row_as_view(MatrixType&& matrix, const uint row) {
    return blaze::row(std::forward<MatrixType>(matrix), row);
}

template <typename MatrixType>
double determinant(MatrixType& matrix) {
    return blaze::det(matrix);
}

template <typename MatrixType>
decltype(auto) inverse(MatrixType& matrix) {
    return blaze::inv(matrix);
}

template <typename T, int m, int n = m, bool SO = blaze::rowMajor>
StatVector<T, m * n> flatten(const StatMatrix<T, m, n>& matrix) {
    StatVector<T, m * n> ret(matrix.data());
    blaze::CustomMatrix<T, blaze::unaligned, blaze::unpadded, SO>(ret.data(), m, n) = matrix;
    return ret;
}

template <typename T, int m, bool SO = blaze::rowMajor>
DynVector<T> flatten(const HybMatrix<T, m>& matrix) {
    uint n = blaze::columns(matrix);

    DynVector<T> ret(m * n, matrix.data());
    blaze::CustomMatrix<T, blaze::unaligned, blaze::unpadded, SO>(ret.data(), m, n) = matrix;
    return ret;
}

template <typename T, bool SO = blaze::rowMajor>
DynVector<T> flatten(const DynMatrix<T>& matrix) {
    uint m = blaze::rows(matrix);
    uint n = blaze::columns(matrix);

    DynVector<T> ret(m * n);
    blaze::CustomMatrix<T, blaze::unaligned, blaze::unpadded, SO>(ret.data(), m, n) = matrix;
    return ret;
}

template <typename VT, bool SO>
decltype(auto) reverse_rows(const blaze::Matrix<VT,SO>& matrix) {
    size_t n_rows = blaze::rows((~matrix));
    std::vector<size_t> reversed_indices(n_rows);
    for ( uint i = 0; i < n_rows; ++i ) {
        reversed_indices[i] = n_rows - 1 - i;
    }
    return blaze::rows((~matrix),reversed_indices);
}

template <typename VT, bool SO>
decltype(auto) reverse_columns(const blaze::Matrix<VT,SO>& matrix) {
    size_t n_rows = blaze::columns((~matrix));
    std::vector<size_t> reversed_indices(n_rows);
    for ( uint i = 0; i < n_rows; ++i ) {
        reversed_indices[i] = n_rows - 1 - i;
    }
    return blaze::columns((~matrix),reversed_indices);
}


template <typename LeftMatrixType, typename RightMatrixType>
decltype(auto) mat_cw_mult(const LeftMatrixType& matrix_left, const RightMatrixType& matrix_right) {
    return matrix_left % matrix_right;
}


template <typename LeftMatrixType, typename RightMatrixType>
decltype(auto) mat_cw_div(const LeftMatrixType& matrix_left, const RightMatrixType& matrix_right) {
    return blaze::map( matrix_left, matrix_right, blaze::Div{} );
}


/* Solving Linear System */
template <typename MatrixType, typename ArrayType>
void solve_sle(MatrixType& A, ArrayType& B) {
    int ipiv[blaze::columns(A)];
    blaze::gesv(A, B, ipiv);
}

template <typename ArrayType, typename T>
void solve_sle(SparseMatrix<T>& A_sparse, ArrayType& B) {
    // Avoid using this function, use a library with sparse solvers, e.g. Eigen
    // Transforming sparse to dense causes ill conditioned problems
    // Solutions generated here can be rubbish
    printf("No sparse solver in Blaze! Consult use_blaze.hpp!\n");
    abort();

    DynMatrix<double> A_dense = A_sparse;
    int ipiv[blaze::columns(A_dense)];
    blaze::gesv(A_dense, B, ipiv);
}

#endif