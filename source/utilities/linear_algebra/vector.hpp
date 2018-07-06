#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "general_definitions.hpp"

#include <eigen3/Eigen/Dense>

namespace Utilities {
namespace LinearAlgebra {
template <typename Type>
struct Matrix;

template <typename Type>
struct Vector {
  private:
    using VectorType = Eigen::Matrix<Type, Eigen::Dynamic, 1>;

    VectorType vector;

    friend class Matrix<Type>;

  public:
    Vector() = default;
    Vector(const uint n_rows);

    void resize(const uint n_rows);

    uint size();
    uint rows();
    uint columns();

    decltype(auto) transpose();
    decltype(auto) dot(const Vector<Type>& vec);

    /* Overloading Operators */

    Type& operator()(const uint i);
    Type operator()(const uint i) const;

    Vector<Type> operator-();

    Vector<Type>& operator+=(const Vector<Type>& vec);  // vec += vec
    Vector<Type>& operator-=(const Vector<Type>& vec);  // vec -= vec
    Vector<Type>& operator*=(const Type& scalar);       // vec *= scalar
    Vector<Type>& operator/=(const Type& scalar);       // vec /= scalar

    decltype(auto) operator+(const Vector<Type>& vec);  // vec + vec
    decltype(auto) operator-(const Vector<Type>& vec);  // vec - vec
    decltype(auto) operator*(const Type& scalar);       // scalar * vec
    decltype(auto) operator/(const Type& scalar);       // vec / scalar

    template <typename T>
    friend decltype(auto) operator*(const T& scalar, const Vector<T>& vec);  // vec * scalar NB: try to avoid

    /* Expression template compatibility */

    template <typename ExpressionType>
    Vector<Type>& operator=(const ExpressionType& expr);

    template <typename ExpressionType>
    Vector<Type>& operator+=(const ExpressionType& expr);  // vec += expr
    template <typename ExpressionType>
    Vector<Type>& operator-=(const ExpressionType& expr);  // vec -= expr
    template <typename ExpressionType>
    Vector<Type>& operator*=(const ExpressionType& expr);  // vec *= expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const ExpressionType& expr, const Vector<T>& vec);  // expr + vec
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const Vector<T>& vec, const ExpressionType& expr);  // vec + expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const ExpressionType& expr, const Vector<T>& vec);  // expr - vec
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const Vector<T>& vec, const ExpressionType& expr);  // vec - expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const Vector<T>& vec, const ExpressionType& expr);  // vec * expr
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const ExpressionType& expr, const Vector<T>& vec);  // expr * vec
};

template <typename Type>
Vector<Type>::Vector(const uint n_rows) : vector(n_rows) {}

template <typename Type>
void Vector<Type>::resize(const uint n_rows) {
    this->vector.conservativeResize(n_rows);
}

template <typename Type>
uint Vector<Type>::size() {
    return this->vector.size();
}

template <typename Type>
uint Vector<Type>::rows() {
    return this->vector.rows();
}

template <typename Type>
uint Vector<Type>::columns() {
    return this->vector.cols();
}

template <typename Type>
decltype(auto) Vector<Type>::transpose() {
    return this->vector.transpose();
}

template <typename Type>
decltype(auto) Vector<Type>::dot(const Vector<Type>& vec) {
    return this->vector.dot(vec.vector);
}

template <typename Type>
Type& Vector<Type>::operator()(const uint i) {
    return this->vector(i);
}

template <typename Type>
Type Vector<Type>::operator()(const uint i) const {
    return this->vector(i);
}

template <typename Type>
Vector<Type> Vector<Type>::operator-() {
    Vector<Type> vec;
    vec.vector = -this->vector;
    return vec;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator+=(const Vector<Type>& vec) {
    this->vector += vec.vector;
    return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator-=(const Vector<Type>& vec) {
    this->vector -= vec.vector;
    return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator*=(const Type& scalar) {
    this->vector *= scalar;
    return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator/=(const Type& scalar) {
    this->vector *= scalar;
    return *this;
}

template <typename Type>
decltype(auto) Vector<Type>::operator+(const Vector<Type>& vec) {
    return this->vector + vec.vector;
}

template <typename Type>
decltype(auto) Vector<Type>::operator-(const Vector<Type>& vec) {
    return this->vector - vec.vector;
}

template <typename Type>
decltype(auto) Vector<Type>::operator*(const Type& scalar) {
    return this->vector * scalar;
}

template <typename Type>
decltype(auto) Vector<Type>::operator/(const Type& scalar) {
    return this->vector / scalar;
}

template <typename T>
decltype(auto) operator*(const T& scalar, const Vector<T>& vec) {
    return scalar * vec.vector;
}

template <typename Type>
template <typename ExpressionType>
Vector<Type>& Vector<Type>::operator=(const ExpressionType& expr) {
    this->vector = expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
Vector<Type>& Vector<Type>::operator+=(const ExpressionType& expr) {
    this->vector += expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
Vector<Type>& Vector<Type>::operator-=(const ExpressionType& expr) {
    this->vector -= expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
Vector<Type>& Vector<Type>::operator*=(const ExpressionType& expr) {
    this->vector *= expr;
    return *this;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator+(const ExpressionType& expr, const Vector<Type>& vec) {
    return expr + vec.vector;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator+(const Vector<Type>& vec, const ExpressionType& expr) {
    return vec.vector + expr;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator-(const ExpressionType& expr, const Vector<Type>& vec) {
    return expr - vec.vector;
}

template <typename ExpressionType, typename Type>
decltype(auto) operator-(const Vector<Type>& vec, const ExpressionType& expr) {
    return vec.vector - expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator*(const Vector<T>& vec, const ExpressionType& expr) {
    return vec.vector * expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator*(const ExpressionType& expr, const Vector<T>& vec) {
    return expr * vec.vector;
}
}
}

#endif