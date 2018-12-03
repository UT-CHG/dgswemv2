#ifndef IS_SOA_HPP
#define IS_SOA_HPP

#include <type_traits>

namespace Utilities {

namespace detail {
template <typename T, typename = void>
struct has_AccessorType : std::false_type {};

template <typename T>
struct has_AccessorType<
    T,
    typename std::enable_if<std::is_object<typename T::AccessorType>::value>::type>
    : std::true_type {};
}


template <typename T, typename = void>
struct is_SoA : std::false_type {};

template <typename T>
struct is_SoA<
    T,
    typename std::enable_if<detail::has_AccessorType<T>::value &&
                            std::is_same<typename T::AccessorType,
                                         decltype(std::declval<T>().at(0u))>::value
                            >::type
    > : std::true_type {};
}

#endif