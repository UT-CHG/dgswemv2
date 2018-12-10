#ifndef IS_VECTORIZED_HPP
#define IS_VECTORIZED_HPP

#include <type_traits>

namespace Utilities {

template <typename T, typename = int>
struct is_vectorized : std::false_type {};

template <typename T>
struct is_vectorized<
    T,
    decltype( (void) T::is_vectorized, 0)
    >
    : std::integral_constant<bool, T::is_vectorized> {};

}

#endif