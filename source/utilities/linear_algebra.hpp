#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <utility>

#ifdef USE_BLAZE
#include "utilities/linear_algebra/use_blaze.hpp"
#endif

#ifdef USE_EIGEN
#include "utilities/linear_algebra/use_eigen.hpp"
#endif

template <>
struct Result<double> {
    using type = double;
};

// The following are STL containers with aligned allocators.
// These should be used whenever the template parameter is
// a static or Hybrid vector type or contains is a class
// which contains Hybrid or static vector types
template <typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;

// On Macbook one can get error:
// static assertion failed: std::map must have the same value_type as its allocator
// static_assert(is_same<typename _Alloc::value_type, value_type>::value
// A fix is adding const to Key in std::pair in AllignedAllocator
// https://github.com/JakobEngel/dso/issues/111
template <typename Key, typename T, typename Compare = std::less<Key>>
using AlignedMap = std::map<Key, T, Compare, AlignedAllocator<std::pair<const Key, T>>>;

namespace detail {
template<typename T, bool SO, size_t... Is>
std::array<DynView<T,SO>,sizeof...(Is)> make_rows_as_views_impl( std::array<DynMatrix<T,SO>,sizeof...(Is)>& mats, size_t index, std::index_sequence<Is...>) {
    return { row_as_view(mats[Is],index)... };
}
}
template<typename T, bool SO, size_t N>
std::array<DynView<T,SO>,N> make_rows_as_views( std::array<DynMatrix<T,SO>,N>& mats, size_t index) {
    using Indices = std::make_index_sequence<N>;
    return detail::make_rows_as_views_impl(mats, index, Indices{});
}
#endif