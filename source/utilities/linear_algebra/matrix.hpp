#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "general_definitions.hpp"
#include <eigen3/Eigen/Dense>

#include "vector.hpp"

namespace Utilities {
namespace LinearAlgebra {
template <typename Type, uint n_rows>
struct Vector;

template <typename Type, uint n_rows, uint n_cols>
struct Matrix {
  private:
    using MatrixType = Eigen::Matrix<Type, n_rows, n_cols>;

    MatrixType matrix;

  public:
    Matrix() = default;

    uint size();
    uint rows();
    uint columns();

    decltype(auto) transpose();

    /* Mat-Mat, Mat-Vec, Vec-Mat(if Mat is 1 row giving outer product) */
    template <uint n_rows_rhs>
    decltype(auto) operator*(const Matrix<Type, n_cols, n_rows_rhs>& mat);  // mat * mat

    decltype(auto) operator*(const Vector<Type, n_cols>& vec);  // mat * vec

    template <typename T, uint nr>
    friend decltype(auto) operator*(const Vector<T, nr>& vec, const Matrix<T, 1, nr>& mat);  // vec * mat(1 row)

    /* Overloading Operators */

    Type& operator()(const uint i, const uint j);
    Type operator()(const uint i, const uint j) const;

    Matrix<Type, n_rows, n_cols> operator-();

    Matrix<Type, n_rows, n_cols>& operator+=(const Matrix<Type, n_rows, n_cols>& mat);  // mat += mat
    Matrix<Type, n_rows, n_cols>& operator-=(const Matrix<Type, n_rows, n_cols>& mat);  // mat -= mat
    Matrix<Type, n_rows, n_cols>& operator*=(const Type& scalar);                       // mat *= scalar
    Matrix<Type, n_rows, n_cols>& operator/=(const Type& scalar);                       // mat /= scalar

    decltype(auto) operator+(const Matrix<Type, n_rows, n_cols>& mat);  // mat + mat
    decltype(auto) operator-(const Matrix<Type, n_rows, n_cols>& mat);  // mat - mat
    decltype(auto) operator*(const Type& scalar);                       // mat * scalar NB: try to avoid
    decltype(auto) operator/(const Type& scalar);                       // mat / scalar

    template <typename T, uint nr, uint nc>
    friend decltype(auto) operator*(const T& scalar, const Matrix<T, nr, nc>& mat);  // scalar * mat NB: try to avoid

    /* Expression template compatibility */

    template <typename ExpressionType>
    Matrix<Type, n_rows, n_cols>& operator=(const ExpressionType& expr);

    template <typename ExpressionType>
    Matrix<Type, n_rows, n_cols>& operator+=(const ExpressionType& expr);  // mat += expr
    template <typename ExpressionType>
    Matrix<Type, n_rows, n_cols>& operator-=(const ExpressionType& expr);  // mat -= expr
    template <typename ExpressionType>
    Matrix<Type, n_rows, n_cols>& operator*=(const ExpressionType& expr);  // mat *= expr

    template <typename ExpressionType, typename T, uint nr, uint nc>
    friend decltype(auto) operator+(const ExpressionType& expr, const Matrix<T, nr, nc>& mat);  // expr + mat
    template <typename ExpressionType, typename T, uint nr, uint nc>
    friend decltype(auto) operator+(const Matrix<T, nr, nc>& mat, const ExpressionType& expr);  // mat + expr

    template <typename ExpressionType, typename T, uint nr, uint nc>
    friend decltype(auto) operator-(const ExpressionType& expr, const Matrix<T, nr, nc>& mat);  // expr - mat
    template <typename ExpressionType, typename T, uint nr, uint nc>
    friend decltype(auto) operator-(const Matrix<T, nr, nc>& mat, const ExpressionType& expr);  // mat - expr

    template <typename ExpressionType, typename T, uint nr, uint nc>
    friend decltype(auto) operator*(const ExpressionType& expr, const Matrix<T, nr, nc>& mat);  // expr * mat
    template <typename ExpressionType, typename T, uint nr, uint nc>
    friend decltype(auto) operator*(const Matrix<T, nr, nc>& mat, const ExpressionType& expr);  // mat * expr
};

template <typename Type, uint n_rows, uint n_cols>
uint Matrix<Type, n_rows, n_cols>::size() {
    return this->matrix.size();
}

template <typename Type, uint n_rows, uint n_cols>
uint Matrix<Type, n_rows, n_cols>::rows() {
    return this->matrix.rows();
}

template <typename Type, uint n_rows, uint n_cols>
uint Matrix<Type, n_rows, n_cols>::columns() {
    return this->matrix.cols();
}

template <typename Type, uint n_rows, uint n_cols>
decltype(auto) Matrix<Type, n_rows, n_cols>::transpose() {
    return this->matrix.transpose();
}

template <typename Type, uint n_rows, uint n_cols>
template <uint n_rows_rhs>
decltype(auto) Matrix<Type, n_rows, n_cols>::operator*(const Matrix<Type, n_cols, n_rows_rhs>& mat) {
    return this->matrix * mat.matrix;
}

template <typename Type, uint n_rows, uint n_cols>
decltype(auto) Matrix<Type, n_rows, n_cols>::operator*(const Vector<Type, n_cols>& vec) {
    return this->matrix * vec.vector;
}

template <typename T, uint nr>
decltype(auto) operator*(const Vector<T, nr>& vec, const Matrix<T, 1, nr>& mat) {
    return vec.vector * mat.matrix;
}

template <typename Type, uint n_rows, uint n_cols>
Type& Matrix<Type, n_rows, n_cols>::operator()(const uint i, const uint j) {
    return this->matrix(i, j);
}

template <typename Type, uint n_rows, uint n_cols>
Type Matrix<Type, n_rows, n_cols>::operator()(const uint i, const uint j) const {
    return this->matrix(i, j);
}

template <typename Type, uint n_rows, uint n_cols>
Matrix<Type, n_rows, n_cols> Matrix<Type, n_rows, n_cols>::operator-() {
    Matrix<Type, n_rows, n_cols> mat;
    mat.matrix = -this->matrix;
    return mat;
}

template <typename Type, uint n_rows, uint n_cols>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator+=(const Matrix<Type, n_rows, n_cols>& mat) {
    this->matrix += mat.matrix;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator-=(const Matrix<Type, n_rows, n_cols>& mat) {
    this->matrix -= mat.matrix;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator*=(const Type& scalar) {
    this->matrix *= scalar;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator/=(const Type& scalar) {
    this->matrix /= scalar;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
decltype(auto) Matrix<Type, n_rows, n_cols>::operator+(const Matrix<Type, n_rows, n_cols>& mat) {
    return this->matrix + mat.matrix;
}

template <typename Type, uint n_rows, uint n_cols>
decltype(auto) Matrix<Type, n_rows, n_cols>::operator-(const Matrix<Type, n_rows, n_cols>& mat) {
    return this->matrix - mat.matrix;
}

template <typename Type, uint n_rows, uint n_cols>
decltype(auto) Matrix<Type, n_rows, n_cols>::operator*(const Type& scalar) {
    return this->matrix * scalar;
}

template <typename Type, uint n_rows, uint n_cols>
decltype(auto) Matrix<Type, n_rows, n_cols>::operator/(const Type& scalar) {
    return this->matrix / scalar;
}

template <typename T, uint nr, uint nc>
decltype(auto) operator*(const T& scalar, const Matrix<T, nr, nc>& mat) {
    return scalar * mat.matrix;
}

template <typename Type, uint n_rows, uint n_cols>
template <typename ExpressionType>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator=(const ExpressionType& expr) {
    this->matrix = expr;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
template <typename ExpressionType>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator+=(const ExpressionType& expr) {
    this->matrix += expr;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
template <typename ExpressionType>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator-=(const ExpressionType& expr) {
    this->matrix -= expr;
    return *this;
}

template <typename Type, uint n_rows, uint n_cols>
template <typename ExpressionType>
Matrix<Type, n_rows, n_cols>& Matrix<Type, n_rows, n_cols>::operator*=(const ExpressionType& expr) {
    this->matrix *= expr;
    return *this;
}

template <typename ExpressionType, typename T, uint nr, uint nc>
decltype(auto) operator+(const ExpressionType& expr, const Matrix<T, nr, nc>& mat) {
    return expr + mat.matrix;
}

template <typename ExpressionType, typename T, uint nr, uint nc>
decltype(auto) operator+(const Matrix<T, nr, nc>& mat, const ExpressionType& expr) {
    return mat.matrix + expr;
}

template <typename ExpressionType, typename T, uint nr, uint nc>
decltype(auto) operator-(const ExpressionType& expr, const Matrix<T, nr, nc>& mat) {
    return expr - mat.matrix;
}

template <typename ExpressionType, typename T, uint nr, uint nc>
decltype(auto) operator-(const Matrix<T, nr, nc>& mat, const ExpressionType& expr) {
    return mat.matrix - expr;
}

template <typename ExpressionType, typename T, uint nr, uint nc>
decltype(auto) operator*(const ExpressionType& expr, const Matrix<T, nr, nc>& mat) {
    return expr * mat.matrix;
}

template <typename ExpressionType, typename T, uint nr, uint nc>
decltype(auto) operator*(const Matrix<T, nr, nc>& mat, const ExpressionType& expr) {
    return mat.matrix * expr;
}
}
}

#endif