#ifndef HETEROGENEOUS_CONTAINERS_HPP
#define HETEROGENEOUS_CONTAINERS_HPP

#include "../general_definitions.hpp"
#include "tuple_helpers.hpp"

namespace Utilities {
/**
 * Vector-like container for storing heterogeneous classes.
 * The implementation of the data is std::tuple<std::vector<Ts>...>
 *
 * @tparam Ts... Types to be stored in the vector
 * @note HeterogeneousVector is not a random access container. Due to retention of strong
 *       typing, access to elements must be accompanied by the desired type.
 */
template <typename... Ts>
struct HeterogeneousVector {
    using TupleType = std::tuple<Ts...>;
    std::tuple<std::vector<Ts>...> data;

    /**
     * Returns the total number of elements in the HeterogeneousVector
     */
    uint size() {
        uint size = 0;

        for_each_in_tuple(this->data, [&size](const auto& vector) { size += vector.size(); });

        return size;
    }

    /**
     * Constructs T(args...) at the end of corresponding vector
     *
     * @tparam T type to be emplaced back
     * @param args... Arguments supplied to constructor of T
     */
    template <typename T, typename... Args>
    void emplace_back(Args... args) {
        static_assert(has_type<T, TupleType>::value, "Error in HeterogeneousVector::emplace_back: Type not found");

        std::get<index<T, TupleType>::value>(data).emplace_back(std::forward<Args>(args)...);
    }

    /**
     * Returns the i-th element of type T with bounds checking
     *
     * @tparam T type of entry to be returned
     * @param i index of the vector of type T to be returned
     */
    template <typename T>
    T& at(uint i) {
        static_assert(has_type<T, TupleType>::value, "Error in HeterogeneousVector::at: Type not found");

        return std::get<index<T, TupleType>::value>(data).at(i);
    }

    /**
     * const at implementation
     */
    template <typename T>
    const T& at(uint i) const {
        static_assert(has_type<T, TupleType>::value, "Error in HeterogeneousVector::at: Type not found");

        return std::get<index<T, TupleType>::value>(data).at(i);
    }
};

/**
 * map-like container for storing heterogeneous classes with unsigned int key.
 * The implementation of the data is std::tuple<std::map<uint,Ts>...>
 *
 * @tparam Ts... Types to be stored in the vector
 */
template <typename... Ts>
struct HeterogeneousMap {
    using TupleType = std::tuple<Ts...>;
    std::tuple<std::map<uint, Ts>...> data;

    /**
     * Returns the total number of elements in the HeterogeneousMap
     */
    uint size() {
        uint size = 0;

        for_each_in_tuple(this->data, [&size](const auto& map) { size += map.size(); });

        return size;
    }

    /**
     * Constructs T(args...) in the corresponding map
     *
     * @tparam T type to be emplaced back
     * @param args... Arguments supplied to constructor of T
     */
    template <typename T>
    void emplace(uint n, T&& t) {
        static_assert(has_type<T, TupleType>::value, "Error in HeterogeneousMap::emplace: Type not found");

        std::get<index<T, TupleType>::value>(data).emplace(std::make_pair(n, std::forward<T>(t)));
    }

    /**
     * Returns the value assocaited with the map of type with bounds checking T
     *
     * @tparam T type of entry to be returned
     * @param key key-value for map
     */
    template <typename T>
    T& at(uint key) {
        static_assert(has_type<T, TupleType>::value, "Error in HeterogeneousMap::at: Type not found");

        return std::get<index<T, TupleType>::value>(data).at(key);
    }

    /**
     * const at implementation
     */
    template <typename T>
    const T& at(uint key) const {
        static_assert(has_type<T, TupleType>::value, "Error in HeterogeneousMap::at: Type not found");

        return std::get<index<T, TupleType>::value>(data).at(key);
    }
};
}

#endif