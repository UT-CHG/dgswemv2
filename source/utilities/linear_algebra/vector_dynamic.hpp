#ifndef VECTOR_DYNAMIC_HPP
#define VECTOR_DYNAMIC_HPP

#include "general_definitions.hpp"

#include <eigen3/Eigen/Dense>

namespace Utilities {
namespace LinearAlgebra {
template <typename Type>
struct MatrixDyn;

template <typename Type>
struct VectorDyn {
  private:
    using VectorType = Eigen::Matrix<Type, Eigen::Dynamic, 1>;

    VectorType vector;

    template <typename T>
    friend struct MatrixDyn;

  public:
    VectorDyn() = default;
    VectorDyn(const uint n_rows);

    void resize(const uint n_rows);

    uint size();
    uint rows();
    uint columns();

    decltype(auto) transpose();
    decltype(auto) dot(const VectorDyn<Type>& vec);

    /* redeclared here defined in matrix */
    template <typename T>
    friend decltype(auto) operator*(const VectorDyn<T>& vec, const MatrixDyn<T>& mat);  // vec * mat(1 row)

    /* Overloading Operators */

    Type& operator()(const uint i);
    Type operator()(const uint i) const;

    VectorDyn<Type> operator-();

    VectorDyn<Type>& operator+=(const VectorDyn<Type>& vec);  // vec += vec
    VectorDyn<Type>& operator-=(const VectorDyn<Type>& vec);  // vec -= vec
    VectorDyn<Type>& operator*=(const Type& scalar);          // vec *= scalar
    VectorDyn<Type>& operator/=(const Type& scalar);          // vec /= scalar

    decltype(auto) operator+(const VectorDyn<Type>& vec);  // vec + vec
    decltype(auto) operator-(const VectorDyn<Type>& vec);  // vec - vec
    decltype(auto) operator*(const Type& scalar);          // scalar * vec
    decltype(auto) operator/(const Type& scalar);          // vec / scalar

    template <typename T>
    friend decltype(auto) operator*(const T& scalar, const VectorDyn<T>& vec);  // vec * scalar NB: try to avoid

    /* Expression template compatibility */

    template <typename ExpressionType>
    VectorDyn<Type>& operator=(const ExpressionType& expr);

    template <typename ExpressionType>
    VectorDyn<Type>& operator+=(const ExpressionType& expr);  // vec += expr
    template <typename ExpressionType>
    VectorDyn<Type>& operator-=(const ExpressionType& expr);  // vec -= expr
    template <typename ExpressionType>
    VectorDyn<Type>& operator*=(const ExpressionType& expr);  // vec *= expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const ExpressionType& expr, const VectorDyn<T>& vec);  // expr + vec
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator+(const VectorDyn<T>& vec, const ExpressionType& expr);  // vec + expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const ExpressionType& expr, const VectorDyn<T>& vec);  // expr - vec
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator-(const VectorDyn<T>& vec, const ExpressionType& expr);  // vec - expr

    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const VectorDyn<T>& vec, const ExpressionType& expr);  // vec * expr
    template <typename ExpressionType, typename T>
    friend decltype(auto) operator*(const ExpressionType& expr, const VectorDyn<T>& vec);  // expr * vec
};

template <typename Type>
VectorDyn<Type>::VectorDyn(const uint n_rows) : vector(n_rows) {}

template <typename Type>
void VectorDyn<Type>::resize(const uint n_rows) {
    this->vector.conservativeResize(n_rows);
}

template <typename Type>
uint VectorDyn<Type>::size() {
    return this->vector.size();
}

template <typename Type>
uint VectorDyn<Type>::rows() {
    return this->vector.rows();
}

template <typename Type>
uint VectorDyn<Type>::columns() {
    return this->vector.cols();
}

template <typename Type>
decltype(auto) VectorDyn<Type>::transpose() {
    return this->vector.transpose();
}

template <typename Type>
decltype(auto) VectorDyn<Type>::dot(const VectorDyn<Type>& vec) {
    return this->vector.dot(vec.vector);
}

template <typename Type>
Type& VectorDyn<Type>::operator()(const uint i) {
    return this->vector(i);
}

template <typename Type>
Type VectorDyn<Type>::operator()(const uint i) const {
    return this->vector(i);
}

template <typename Type>
VectorDyn<Type> VectorDyn<Type>::operator-() {
    VectorDyn<Type> vec;
    vec.vector = -this->vector;
    return vec;
}

template <typename Type>
VectorDyn<Type>& VectorDyn<Type>::operator+=(const VectorDyn<Type>& vec) {
    this->vector += vec.vector;
    return *this;
}

template <typename Type>
VectorDyn<Type>& VectorDyn<Type>::operator-=(const VectorDyn<Type>& vec) {
    this->vector -= vec.vector;
    return *this;
}

template <typename Type>
VectorDyn<Type>& VectorDyn<Type>::operator*=(const Type& scalar) {
    this->vector *= scalar;
    return *this;
}

template <typename Type>
VectorDyn<Type>& VectorDyn<Type>::operator/=(const Type& scalar) {
    this->vector *= scalar;
    return *this;
}

template <typename Type>
decltype(auto) VectorDyn<Type>::operator+(const VectorDyn<Type>& vec) {
    return this->vector + vec.vector;
}

template <typename Type>
decltype(auto) VectorDyn<Type>::operator-(const VectorDyn<Type>& vec) {
    return this->vector - vec.vector;
}

template <typename Type>
decltype(auto) VectorDyn<Type>::operator*(const Type& scalar) {
    return this->vector * scalar;
}

template <typename Type>
decltype(auto) VectorDyn<Type>::operator/(const Type& scalar) {
    return this->vector / scalar;
}

template <typename T>
decltype(auto) operator*(const T& scalar, const VectorDyn<T>& vec) {
    return scalar * vec.vector;
}

template <typename Type>
template <typename ExpressionType>
VectorDyn<Type>& VectorDyn<Type>::operator=(const ExpressionType& expr) {
    this->vector = expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
VectorDyn<Type>& VectorDyn<Type>::operator+=(const ExpressionType& expr) {
    this->vector += expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
VectorDyn<Type>& VectorDyn<Type>::operator-=(const ExpressionType& expr) {
    this->vector -= expr;
    return *this;
}

template <typename Type>
template <typename ExpressionType>
VectorDyn<Type>& VectorDyn<Type>::operator*=(const ExpressionType& expr) {
    this->vector *= expr;
    return *this;
}

template <typename ExpressionType, typename T>
decltype(auto) operator+(const ExpressionType& expr, const VectorDyn<T>& vec) {
    return expr + vec.vector;
}

template <typename ExpressionType, typename T>
decltype(auto) operator+(const VectorDyn<T>& vec, const ExpressionType& expr) {
    return vec.vector + expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator-(const ExpressionType& expr, const VectorDyn<T>& vec) {
    return expr - vec.vector;
}

template <typename ExpressionType, typename T>
decltype(auto) operator-(const VectorDyn<T>& vec, const ExpressionType& expr) {
    return vec.vector - expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator*(const VectorDyn<T>& vec, const ExpressionType& expr) {
    return vec.vector * expr;
}

template <typename ExpressionType, typename T>
decltype(auto) operator*(const ExpressionType& expr, const VectorDyn<T>& vec) {
    return expr * vec.vector;
}
}
}

#endif