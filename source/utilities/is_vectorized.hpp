#ifndef IS_VECTORIZED_HPP
#define IS_VECTORIZED_HPP

#include <type_traits>

namespace Utilities {

template <typename T, typename AccessorType = void, typename = void>
struct is_vectorized : std::false_type {};

template <typename T, typename AccessorType>
struct is_vectorized<
    T, AccessorType,
    typename std::enable_if<
        T::template is_vectorized<AccessorType>()
        >::type
    > : std::integral_constant<bool, T::template is_vectorized<AccessorType>()> {};

}

#endif