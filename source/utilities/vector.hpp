#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "general_definitions.hpp"

#include <eigen3/Eigen/Dense>

namespace Utilities {
namespace LinearAlgebra {
template <typename T>
struct Vector {
  private:
    using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    VectorType vec;

  public:
    Vector() = default;
    Vector(const uint n_rows);

    void resize(const uint n_rows);

    uint size();
    uint rows();
    uint columns();

    T& operator[](const uint i);
    T operator[](const uint i) const;

    template <typename ExpressionType>
    Vector<T>& operator=(const ExpressionType& expr);

    decltype(auto) operator+(const Vector<T>& rvec);
    decltype(auto) operator-(const Vector<T>& rvec);
    decltype(auto) operator*(const T& scalar);
    decltype(auto) operator/(const T& scalar);

    template <typename ExpressionType, typename Type>
    friend decltype(auto) operator+(const ExpressionType& expr, const Vector<Type> rvec);
    template <typename ExpressionType, typename Type>
    friend decltype(auto) operator-(const ExpressionType& expr, const Vector<Type> rvec);
    template <typename Type>
    friend decltype(auto) operator*(const Type& scalar, const Vector<Type> rvec);
    template <typename Type>
    friend decltype(auto) operator/(const Type& scalar, const Vector<Type> rvec);
};

template <typename T>
Vector<T>::Vector(const uint n_rows) : vec(n_rows) {}

template <typename T>
void Vector<T>::resize(const uint n_rows) {
    this->vec.conservativeResize(n_rows);
}

template <typename T>
uint Vector<T>::size() {
    return this->vec.size();
}

template <typename T>
uint Vector<T>::rows() {
    return this->vec.rows();
}

template <typename T>
uint Vector<T>::columns() {
    return this->vec.cols();
}

template <typename T>
T& Vector<T>::operator[](const uint i) {
    return this->vec(i);
}

template <typename T>
T Vector<T>::operator[](const uint i) const {
    return this->vec(i);
}

template <typename T>
template <typename ExpressionType>
Vector<T>& Vector<T>::operator=(const ExpressionType& expr) {
    this->vec = expr;
    return *this;
}

template <typename T>
decltype(auto) Vector<T>::operator+(const Vector<T>& rvec) {
    return this->vec + rvec.vec;
}

template <typename ExpressionType, typename T>
decltype(auto) operator+(const ExpressionType& expr, const Vector<T> rvec) {
    return expr + rvec.vec;
}

template <typename T>
decltype(auto) Vector<T>::operator-(const Vector<T>& rvec) {
    return this->vec - rvec.vec;
}

template <typename ExpressionType, typename T>
decltype(auto) operator-(const ExpressionType& expr, const Vector<T> rvec) {
    return expr - rvec.vec;
}

template <typename T>
decltype(auto) Vector<T>::operator*(const T& scalar) {
    return this->vec * scalar;
}

template <typename T>
decltype(auto) Vector<T>::operator/(const T& scalar) {
    return this->vec / scalar;
}

template <typename Type>
decltype(auto) operator*(const Type& scalar, const Vector<Type> rvec) {
    return rvec.vec * scalar;
}

template <typename Type>
decltype(auto) operator/(const Type& scalar, const Vector<Type> rvec) {
    return rvec.vec / scalar;
}
}
}

#endif