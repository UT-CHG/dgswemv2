#ifndef AT_EACH_HPP
#define AT_EACH_HPP

#include "utilities/linear_algebra.hpp"
#include "utilities/is_soa.hpp"

#include <utility>

namespace Utilities {

namespace detail {
template <typename SoA,
          size_t... Is>
std::array<typename SoA::AccessorType, sizeof...(Is)> at_each_array_impl(std::array<SoA,sizeof...(Is)>& arr,
                                                                         uint index,
                                                                         std::index_sequence<Is...>) {
    return { arr[Is].at(index)... };
}
}

template <typename SoA,
          size_t N,
          typename = typename std::enable_if<is_SoA<SoA>::value>::type
          >
std::array<typename SoA::AccessorType, N> at_each(std::array<SoA,N>& arr, uint index) {
    using Indices = std::make_index_sequence<N>;
    return detail::at_each_array_impl(arr, index, Indices{});
}


template <typename T, template<typename> typename Allocator>
std::vector<std::reference_wrapper<T>, Allocator<std::reference_wrapper<T>>> at_each(std::vector<std::vector<T>,Allocator<std::vector<T>>>& v, uint index) {

    assert( std::all_of( v.cbegin(), v.cend(), [index](const std::vector<T>& sub_v) {
                return index < sub_v.size(); })
        );
    std::vector<std::reference_wrapper<T>, Allocator<std::reference_wrapper<T>>> refs;
    refs.reserve(v.size());
    for ( uint i = 0; i < v.size(); ++i ) {
        refs.push_back( v[i][index] );
    }

    return refs;
}

template <
    typename SoA,
    template<typename> typename Allocator,
    typename = typename std::enable_if<is_SoA<SoA>::value>::type
    >
std::vector<typename SoA::AccessorType, Allocator<typename SoA::AccessorType>> at_each(std::vector<SoA, Allocator<SoA>>& v, uint index) {
    std::vector<typename SoA::AccessorType, Allocator<typename SoA::AccessorType>> refs;
    refs.reserve(v.size());
    for ( uint i = 0; i < v.size(); ++i ) {
        refs.push_back( v[i].at(index) );
    }

    return refs;
}

}


#endif