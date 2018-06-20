#ifndef TUPLE_HELPERS_HPP
#define TUPLE_HELPERS_HPP

#include "general_definitions.hpp"
#include "ignore.hpp"

namespace Utilities {
///////////////////////////////////////////////////////////////////////////////////////////////
// has_type implementation from
// https://stackoverflow.com/questions/25958259/how-do-i-find-out-if-a-tuple-contains-a-type
template <typename T, typename Tuple>
struct has_type;

template <typename T>
struct has_type<T, std::tuple<>> : std::false_type {};

template <typename T, typename U, typename... Ts>
struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

/**
 * Metaprogramming function that statically determines whether a type is in a tuple.
 * has_type inherits value_type based whether or not T is found in Ts...
 *
 * @tparam T the type were looking for
 * @tparam Ts... Parameter pack in which we are looking for T
 */
template <typename T, typename... Ts>
struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

///////////////////////////////////////////////////////////////////////////////////////////////
// index implementation
template <typename T, typename Tuple>
struct index;

/**
 * index returns the first location of a type
 * Based on
 * <a href="https://stackoverflow.com/a/18063608">Casey's implementation</a>.
 *
 * @tparam T type being looked for
 * @tparam Ts... template pack in which we are searching for T
 * @warning There is no current bounds checking for if T is not in Ts... Therefore we recommend
 *          using a static assert with Utilities::has_type to determine whether index is a
 *          well-defined class
 */
template <typename T, typename... Ts>
struct index<T, std::tuple<T, Ts...>> {
    /**
     * Compile time constant of the index of type T in pack Ts...
     */
    static const std::size_t value = 0;
};

template <typename T, typename U, typename... Ts>
struct index<T, std::tuple<U, Ts...>> {
    static const std::size_t value = 1 + index<T, std::tuple<Ts...>>::value;
};

///////////////////////////////////////////////////////////////////////////////////////////////
// for_each_in_tuple implementation
template <int... Is>
struct seq {};

template <int N, int... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

template <int... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

template <typename T, typename F, int... Is>
void for_each(T&& t, F f, seq<Is...>) {
    // This warning is pretty verbose given all the metaprogramming.
    // using compiler macros to supress the warnings in relevant compilers
    auto unused = {(f(std::get<Is>(t)), 0)...};
    ignore(unused);
}

/**
 * Apply the unary function F to each member of the tuple std::tuple<Ts...>.
 * This implementation is based on
 * <a href="https://stackoverflow.com/a/16387374">Andy Prowl's implementation</a>.
 *
 * @param t Tuple with members of type Ts...
 * @param f the unary function.
 */
template <typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...>& t, F f) {
    for_each(t, f, gen_seq<sizeof...(Ts)>());
}

template <typename Tup1, typename Tup2>
struct tuple_join;

/**
 * Metafunction for joining two tuples.
 *
 * @tparam Pack1 first tuple of type std::tuple<Pack1...>
 * @tparam Pack2 second tuple of type std::tuple<Pack2...>
 */
template <typename... Pack1, typename... Pack2>
struct tuple_join<std::tuple<Pack1...>, std::tuple<Pack2...>> {
    /**
     *  joined tuple type std::tuple<Pack1..., Pack2...>
     */
    typedef std::tuple<Pack1..., Pack2...> type;
};
}

#endif
