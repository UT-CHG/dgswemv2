#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "general_definitions.hpp"

#include <eigen3/Eigen/Dense>

namespace Utilities {
namespace LinearAlgebra {
template <typename Type, uint n_rows, uint n_cols>
struct Matrix;

template <typename Type, uint n_rows>
struct Vector {
  private:
    using VectorType = Eigen::Matrix<Type, n_rows, 1>;

    VectorType vector;

    template <typename T, uint nr, uint nc>
    friend struct Matrix;

  public:
    Vector() = default;

    uint size();
    uint rows();
    uint columns();

    decltype(auto) transpose();
    decltype(auto) dot(const Vector<Type, n_rows>& vec);

    /* redeclared here defined in matrix */
    template <typename T, uint nr>
    friend decltype(auto) operator*(const Vector<T, nr>& vec, const Matrix<T, 1, nr>& mat);  // vec * mat(1 row)

    /* Overloading Operators */

    Type& operator()(const uint i);
    Type operator()(const uint i) const;

    Vector<Type, n_rows> operator-();

    Vector<Type, n_rows>& operator+=(const Vector<Type, n_rows>& vec);  // vec += vec
    Vector<Type, n_rows>& operator-=(const Vector<Type, n_rows>& vec);  // vec -= vec
    Vector<Type, n_rows>& operator*=(const Type& scalar);               // vec *= scalar
    Vector<Type, n_rows>& operator/=(const Type& scalar);               // vec /= scalar

    decltype(auto) operator+(const Vector<Type, n_rows>& vec);  // vec + vec
    decltype(auto) operator-(const Vector<Type, n_rows>& vec);  // vec - vec
    decltype(auto) operator*(const Type& scalar);               // vec * scalar NB: try to avoid
    decltype(auto) operator/(const Type& scalar);               // vec / scalar

    template <typename T, uint nr>
    friend decltype(auto) operator*(const T& scalar, const Vector<T, nr>& vec);  // scalar * vec

    /* Expression template compatibility */

    template <typename ExpressionType>
    Vector<Type, n_rows>& operator=(const ExpressionType& expr);

    template <typename ExpressionType>
    Vector<Type, n_rows>& operator+=(const ExpressionType& expr);  // vec += expr
    template <typename ExpressionType>
    Vector<Type, n_rows>& operator-=(const ExpressionType& expr);  // vec -= expr
    template <typename ExpressionType>
    Vector<Type, n_rows>& operator*=(const ExpressionType& expr);  // vec *= expr

    template <typename ExpressionType, typename T, uint nr>
    friend decltype(auto) operator+(const ExpressionType& expr, const Vector<T, nr>& vec);  // expr + vec
    template <typename ExpressionType, typename T, uint nr>
    friend decltype(auto) operator+(const Vector<T, nr>& vec, const ExpressionType& expr);  // vec + expr

    template <typename ExpressionType, typename T, uint nr>
    friend decltype(auto) operator-(const ExpressionType& expr, const Vector<T, nr>& vec);  // expr - vec
    template <typename ExpressionType, typename T, uint nr>
    friend decltype(auto) operator-(const Vector<T, nr>& vec, const ExpressionType& expr);  // vec - expr

    template <typename ExpressionType, typename T, uint nr>
    friend decltype(auto) operator*(const Vector<T, nr>& vec, const ExpressionType& expr);  // vec * expr
    template <typename ExpressionType, typename T, uint nr>
    friend decltype(auto) operator*(const ExpressionType& expr, const Vector<T, nr>& vec);  // expr * vec
};

template <typename Type, uint n_rows>
uint Vector<Type, n_rows>::size() {
    return this->vector.size();
}

template <typename Type, uint n_rows>
uint Vector<Type, n_rows>::rows() {
    return this->vector.rows();
}

template <typename Type, uint n_rows>
uint Vector<Type, n_rows>::columns() {
    return this->vector.cols();
}

template <typename Type, uint n_rows>
decltype(auto) Vector<Type, n_rows>::transpose() {
    return this->vector.transpose();
}

template <typename Type, uint n_rows>
decltype(auto) Vector<Type, n_rows>::dot(const Vector<Type, n_rows>& vec) {
    return this->vector.dot(vec.vector);
}

template <typename Type, uint n_rows>
Type& Vector<Type, n_rows>::operator()(const uint i) {
    return this->vector(i);
}

template <typename Type, uint n_rows>
Type Vector<Type, n_rows>::operator()(const uint i) const {
    return this->vector(i);
}

template <typename Type, uint n_rows>
Vector<Type, n_rows> Vector<Type, n_rows>::operator-() {
    Vector<Type, n_rows> vec;
    vec.vector = -this->vector;
    return vec;
}

template <typename Type, uint n_rows>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator+=(const Vector<Type, n_rows>& vec) {
    this->vector += vec.vector;
    return *this;
}

template <typename Type, uint n_rows>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator-=(const Vector<Type, n_rows>& vec) {
    this->vector -= vec.vector;
    return *this;
}

template <typename Type, uint n_rows>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator*=(const Type& scalar) {
    this->vector *= scalar;
    return *this;
}

template <typename Type, uint n_rows>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator/=(const Type& scalar) {
    this->vector *= scalar;
    return *this;
}

template <typename Type, uint n_rows>
decltype(auto) Vector<Type, n_rows>::operator+(const Vector<Type, n_rows>& vec) {
    return this->vector + vec.vector;
}

template <typename Type, uint n_rows>
decltype(auto) Vector<Type, n_rows>::operator-(const Vector<Type, n_rows>& vec) {
    return this->vector - vec.vector;
}

template <typename Type, uint n_rows>
decltype(auto) Vector<Type, n_rows>::operator*(const Type& scalar) {
    return this->vector * scalar;
}

template <typename Type, uint n_rows>
decltype(auto) Vector<Type, n_rows>::operator/(const Type& scalar) {
    return this->vector / scalar;
}

template <typename T, uint nr>
decltype(auto) operator*(const T& scalar, const Vector<T, nr>& vec) {
    return scalar * vec.vector;
}

template <typename Type, uint n_rows>
template <typename ExpressionType>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator=(const ExpressionType& expr) {
    this->vector = expr;
    return *this;
}

template <typename Type, uint n_rows>
template <typename ExpressionType>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator+=(const ExpressionType& expr) {
    this->vector += expr;
    return *this;
}

template <typename Type, uint n_rows>
template <typename ExpressionType>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator-=(const ExpressionType& expr) {
    this->vector -= expr;
    return *this;
}

template <typename Type, uint n_rows>
template <typename ExpressionType>
Vector<Type, n_rows>& Vector<Type, n_rows>::operator*=(const ExpressionType& expr) {
    this->vector *= expr;
    return *this;
}

template <typename ExpressionType, typename T, uint nr>
decltype(auto) operator+(const ExpressionType& expr, const Vector<T, nr>& vec) {
    return expr + vec.vector;
}

template <typename ExpressionType, typename T, uint nr>
decltype(auto) operator+(const Vector<T, nr>& vec, const ExpressionType& expr) {
    return vec.vector + expr;
}

template <typename ExpressionType, typename T, uint nr>
decltype(auto) operator-(const ExpressionType& expr, const Vector<T, nr>& vec) {
    return expr - vec.vector;
}

template <typename ExpressionType, typename T, uint nr>
decltype(auto) operator-(const Vector<T, nr>& vec, const ExpressionType& expr) {
    return vec.vector - expr;
}

template <typename ExpressionType, typename T, uint nr>
decltype(auto) operator*(const Vector<T, nr>& vec, const ExpressionType& expr) {
    return vec.vector * expr;
}

template <typename ExpressionType, typename T, uint nr>
decltype(auto) operator*(const ExpressionType& expr, const Vector<T, nr>& vec) {
    return expr * vec.vector;
}
}
}

#endif