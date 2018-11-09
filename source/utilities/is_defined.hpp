#ifndef UTILITIES_IS_DEFINED_HPP
#define UTILITIES_IS_DEFINED_HPP

#include <type_traits>

namespace Utilities {

template <typename T, typename = void>
struct is_defined : std::false_type {};

// https://stackoverflow.com/a/39816909
template <typename T>
struct is_defined<
    T,
    typename std::enable_if<std::is_object<T>::value && !std::is_pointer<T>::value && (sizeof(T) > 0)>::type>
    : std::true_type {};
}
#endif