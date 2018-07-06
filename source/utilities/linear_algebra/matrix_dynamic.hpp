#ifndef MATRIX_DYNAMIC_HPP
#define MATRIX_DYNAMIC_HPP

#include "general_definitions.hpp"
#include <eigen3/Eigen/Dense>

#include "vector_dynamic.hpp"

namespace Utilities {
namespace LinearAlgebra {
template <typename Type>
struct VectorDyn;

template <typename Type>
struct MatrixDyn {
  private:
    using MatrixType = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;

    MatrixType matrix;

  public:
    MatrixDyn() = default;
    MatrixDyn(const uint n_rows, const uint n_cols);

    void resize(const uint n_rows, const uint n_cols);

    uint size();
    uint rows();
    uint columns();

    decltype(auto) transpose();

    /* Mat-Mat, Mat-Vec, Vec-Mat(if Mat is 1 row giving outer product) */
    decltype(auto) operator*(const MatrixDyn<Type>& mat);  // mat * mat

    decltype(auto) operator*(const VectorDyn<Type>& vec);  // mat * vec

    template <typename T>
    friend decltype(auto) operator*(const VectorDyn<T>& vec, const MatrixDyn<T>& mat);  // vec * mat(1 row)

    /* Overloading Operators */

    Type& operator()(const uint i, const uint j);
    Type operator()(const uint i, const uint j) const;

    MatrixDyn<Type> operator-();

    MatrixDyn<Type>& operator+=(const MatrixDyn<Type>& mat);  // mat += mat
    MatrixDyn<Type>& operator-=(const MatrixDyn<Type>& mat);  // mat -= mat
    MatrixDyn<Type>& operator*=(const Type& scalar);          // mat *= scalar
    MatrixDyn<Type>& operator/=(const Type& scalar);          // mat /= scalar

    decltype(auto) operator+(const MatrixDyn<Type>& mat);  // mat + mat
    decltype(auto) operator-(const MatrixDyn<Type>& mat);  // mat - mat
    decltype(auto) operator*(const Type& scalar);          // scalar * mat
    decltype(auto) operator/(const Type& scalar);          // mat / scalar

    template <typename T>
    friend decltype(auto) operator*(const T& scalar, const MatrixDyn<T>& mat);  // mat * scalar NB: try to avoid

    /* Expression template compatibility */

    template <typename ExpressionType>
    MatrixDyn<Type>& operator=(const ExpressionType& expr);

    template <typename ExpressionType>
    MatrixDyn<Type>& operator+=(const ExpressionType& expr);  // mat += expr
    template <typename ExpressionType>
    MatrixDyn<Type>& operator-=(const ExpressionType& expr);  // mat -= expr
    template <typename ExpressionType>
    MatrixDyn<Type>& operator*=(const ExpressionType& expr);  // mat *= expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const ExpressionType& expr, const MatrixDyn<T>& mat);  // expr + mat
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const MatrixDyn<T>& mat, const ExpressionType& expr);  // mat + expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const ExpressionType& expr, const MatrixDyn<T>& mat);  // expr - mat
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const MatrixDyn<T>& mat, const ExpressionType& expr);  // mat - expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const ExpressionType& expr, const MatrixDyn<T>& mat);  // expr * mat
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const MatrixDyn<T>& mat, const ExpressionType& expr);  // mat * expr
};

template <typename Type>
MatrixDyn<Type>::MatrixDyn(const uint n_rows, const uint n_cols) : matrix(n_rows, n_cols) {}

template <typename Type>
void MatrixDyn<Type>::resize(const uint n_rows, const uint n_cols) {
    this->matrix.conservativeResize(n_rows, n_cols);
}

template <typename Type>
uint MatrixDyn<Type>::size() {
    return this->matrix.size();
}

template <typename Type>
uint MatrixDyn<Type>::rows() {
    return this->matrix.rows();
}

template <typename Type>
uint MatrixDyn<Type>::columns() {
    return this->matrix.cols();
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::transpose() {
    return this->matrix.transpose();
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::operator*(const MatrixDyn<Type>& mat) {
    return this->matrix * mat.matrix;
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::operator*(const VectorDyn<Type>& vec) {
    return this->matrix * vec.vector;
}

template <typename T>
decltype(auto) operator*(const VectorDyn<T>& vec, const MatrixDyn<T>& mat) {
    return vec.vector * mat.matrix;
}

template <typename Type>
Type& MatrixDyn<Type>::operator()(const uint i, const uint j) {
    return this->matrix(i, j);
}

template <typename Type>
Type MatrixDyn<Type>::operator()(const uint i, const uint j) const {
    return this->matrix(i, j);
}

template <typename Type>
MatrixDyn<Type> MatrixDyn<Type>::operator-() {
    MatrixDyn<Type> mat;
    mat.matrix = -this->matrix;
    return mat;
}

template <typename Type>
MatrixDyn<Type>& MatrixDyn<Type>::operator+=(const MatrixDyn<Type>& mat) {
    this->matrix += mat.matrix;
    return *this;
}

template <typename Type>
MatrixDyn<Type>& MatrixDyn<Type>::operator-=(const MatrixDyn<Type>& mat) {
    this->matrix -= mat.matrix;
    return *this;
}

template <typename Type>
MatrixDyn<Type>& MatrixDyn<Type>::operator*=(const Type& scalar) {
    this->matrix *= scalar;
    return *this;
}

template <typename Type>
MatrixDyn<Type>& MatrixDyn<Type>::operator/=(const Type& scalar) {
    this->matrix /= scalar;
    return *this;
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::operator+(const MatrixDyn<Type>& mat) {
    return this->matrix + mat.matrix;
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::operator-(const MatrixDyn<Type>& mat) {
    return this->matrix - mat.matrix;
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::operator*(const Type& scalar) {
    return this->matrix * scalar;
}

template <typename Type>
decltype(auto) MatrixDyn<Type>::operator/(const Type& scalar) {
    return this->matrix / scalar;
}

template <typename T>
decltype(auto) operator*(const T& scalar, const MatrixDyn<T>& mat) {
    return scalar * mat.matrix;
}

template <typename Type>
template <typename ExpressionType>
MatrixDyn<Type>& MatrixDyn<Type>::operator=(const ExpressionType& expr) {
    this->matrix = expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
MatrixDyn<Type>& MatrixDyn<Type>::operator+=(const ExpressionType& expr) {
    this->matrix += expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
MatrixDyn<Type>& MatrixDyn<Type>::operator-=(const ExpressionType& expr) {
    this->matrix -= expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
MatrixDyn<Type>& MatrixDyn<Type>::operator*=(const ExpressionType& expr) {
    this->matrix *= expr;
    return *this;
}

template <typename ExpressionType, typename T>
decltype(auto) operator+(const ExpressionType& expr, const MatrixDyn<T>& mat) {
    return expr + mat.matrix;
}

template <typename ExpressionType, typename T>
decltype(auto) operator+(const MatrixDyn<T>& mat, const ExpressionType& expr) {
    return mat.matrix + expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator-(const ExpressionType& expr, const MatrixDyn<T>& mat) {
    return expr - mat.matrix;
}

template <typename ExpressionType, typename T>
decltype(auto) operator-(const MatrixDyn<T>& mat, const ExpressionType& expr) {
    return mat.matrix - expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator*(const ExpressionType& expr, const MatrixDyn<T>& mat) {
    return expr * mat.matrix;
}

template <typename ExpressionType, typename T>
decltype(auto) operator*(const MatrixDyn<T>& mat, const ExpressionType& expr) {
    return mat.matrix * expr;
}
}
}

#endif