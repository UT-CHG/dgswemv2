#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "general_definitions.hpp"
#include <eigen3/Eigen/Dense>

#include "vector.hpp"

namespace Utilities {
namespace LinearAlgebra {
template <typename Type>
struct Vector;

template <typename Type>
struct Matrix {
  private:
    using MatrixType = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;

    MatrixType matrix;

    friend class Vector<Type>;

  public:
    Matrix() = default;
    Matrix(const uint n_rows, const uint n_cols);

    void resize(const uint n_rows, const uint n_cols);

    uint size();
    uint rows();
    uint columns();

    decltype(auto) transpose();

    /* Mat-Mat and Mat-Vec */
    decltype(auto) operator*(const Matrix<Type>& mat);  // mat * mat
    decltype(auto) operator*(const Vector<Type>& vec);  // mat * vec

    /* Overloading Operators */

    Type& operator()(const uint i, const uint j);
    Type operator()(const uint i, const uint j) const;

    Matrix<Type> operator-();

    Matrix<Type>& operator+=(const Matrix<Type>& mat);  // mat += mat
    Matrix<Type>& operator-=(const Matrix<Type>& mat);  // mat -= mat
    Matrix<Type>& operator*=(const Type& scalar);       // mat *= scalar
    Matrix<Type>& operator/=(const Type& scalar);       // mat /= scalar

    decltype(auto) operator+(const Matrix<Type>& mat);  // mat + mat
    decltype(auto) operator-(const Matrix<Type>& mat);  // mat - mat
    decltype(auto) operator*(const Type& scalar);       // scalar * mat
    decltype(auto) operator/(const Type& scalar);       // mat / scalar

    template <typename T>
    friend decltype(auto) operator*(const T& scalar, const Matrix<T>& mat);  // mat * scalar NB: try to avoid

    /* Expression template compatibility */

    template <typename ExpressionType>
    Matrix<Type>& operator=(const ExpressionType& expr);

    template <typename ExpressionType>
    Matrix<Type>& operator+=(const ExpressionType& expr);  // mat += expr
    template <typename ExpressionType>
    Matrix<Type>& operator-=(const ExpressionType& expr);  // mat -= expr
    template <typename ExpressionType>
    Matrix<Type>& operator*=(const ExpressionType& expr);  // mat *= expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const ExpressionType& expr, const Matrix<T>& mat);  // expr + mat
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const Matrix<T>& mat, const ExpressionType& expr);  // mat + expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const ExpressionType& expr, const Matrix<T>& mat);  // expr - mat
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const Matrix<T>& mat, const ExpressionType& expr);  // mat - expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const ExpressionType& expr, const Matrix<T>& mat);  // expr * mat
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const Matrix<T>& mat, const ExpressionType& expr);  // mat * expr
};

template <typename Type>
Matrix<Type>::Matrix(const uint n_rows, const uint n_cols) : matrix(n_rows, n_cols) {}

template <typename Type>
void Matrix<Type>::resize(const uint n_rows, const uint n_cols) {
    this->matrix.conservativeResize(n_rows, n_cols);
}

template <typename Type>
uint Matrix<Type>::size() {
    return this->matrix.size();
}

template <typename Type>
uint Matrix<Type>::rows() {
    return this->matrix.rows();
}

template <typename Type>
uint Matrix<Type>::columns() {
    return this->matrix.cols();
}

template <typename Type>
decltype(auto) Matrix<Type>::transpose() {
    return this->matrix.transpose();
}

template <typename Type>
decltype(auto) Matrix<Type>::operator*(const Matrix<Type>& mat) {
    return this->matrix * mat.matrix;
}

template <typename Type>
decltype(auto) Matrix<Type>::operator*(const Vector<Type>& vec) {
    return this->matrix * vec.vector;
}

template <typename Type>
Type& Matrix<Type>::operator()(const uint i, const uint j) {
    return this->matrix(i, j);
}

template <typename Type>
Type Matrix<Type>::operator()(const uint i, const uint j) const {
    return this->matrix(i, j);
}

template <typename Type>
Matrix<Type> Matrix<Type>::operator-() {
    Matrix<Type> mat;
    mat.matrix = -this->matrix;
    return mat;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=(const Matrix<Type>& mat) {
    this->matrix += mat.matrix;
    return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=(const Matrix<Type>& mat) {
    this->matrix -= mat.matrix;
    return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=(const Type& scalar) {
    this->matrix *= scalar;
    return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=(const Type& scalar) {
    this->matrix /= scalar;
    return *this;
}

template <typename Type>
decltype(auto) Matrix<Type>::operator+(const Matrix<Type>& mat) {
    return this->matrix + mat.matrix;
}

template <typename Type>
decltype(auto) Matrix<Type>::operator-(const Matrix<Type>& mat) {
    return this->matrix - mat.matrix;
}

template <typename Type>
decltype(auto) Matrix<Type>::operator*(const Type& scalar) {
    return this->matrix * scalar;
}

template <typename Type>
decltype(auto) Matrix<Type>::operator/(const Type& scalar) {
    return this->matrix / scalar;
}

template <typename T>
decltype(auto) operator*(const T& scalar, const Matrix<T>& mat) {
    return scalar * mat.matrix;
}

template <typename Type>
template <typename ExpressionType>
Matrix<Type>& Matrix<Type>::operator=(const ExpressionType& expr) {
    this->matrix = expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
Matrix<Type>& Matrix<Type>::operator+=(const ExpressionType& expr) {
    this->matrix += expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
Matrix<Type>& Matrix<Type>::operator-=(const ExpressionType& expr) {
    this->matrix -= expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
Matrix<Type>& Matrix<Type>::operator*=(const ExpressionType& expr) {
    this->matrix *= expr;
    return *this;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator+(const ExpressionType& expr, const Matrix<Type>& mat) {
    return expr + mat.matrix;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator+(const Matrix<Type>& mat, const ExpressionType& expr) {
    return mat.matrix + expr;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator-(const ExpressionType& expr, const Matrix<Type>& mat) {
    return expr - mat.matrix;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator-(const Matrix<Type>& mat, const ExpressionType& expr) {
    return mat.matrix - expr;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator*(const ExpressionType& expr, const Matrix<Type>& mat) {
    return expr * mat.matrix;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator*(const Matrix<Type>& mat, const ExpressionType& expr) {
    return mat.matrix * expr;
}
}
}

#endif