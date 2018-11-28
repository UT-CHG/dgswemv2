#ifndef AT_EACH_HPP
#define AT_EACH_HPP

#include "utilities/linear_algebra.hpp"
#include "utilities/is_soa.hpp"

namespace Utilities {

template <typename T, typename Allocator>
std::vector<std::reference_wrapper<T>> at_each(std::vector<std::vector<T>,Allocator>& v, uint index) {

    assert( std::all_of( v.cbegin(), v.cend(), [index](const std::vector<T>& sub_v) {
                return index < sub_v.size(); })
        );
    std::vector<std::reference_wrapper<T>> refs;
    refs.reserve(v.size());
    for ( uint i = 0; i < v.size(); ++i ) {
        refs.push_back( v[i][index] );
    }

    return refs;
}

template <
    typename SoA,
    typename = typename std::enable_if<is_SoA<SoA>::value>::type
    >
std::vector<typename SoA::AccessorType> at_each(std::vector<SoA>& v, uint index) {
    std::vector<typename SoA::AccessorType> refs;
    refs.reserve(v.size());
    for ( uint i = 0; i < v.size(); ++i ) {
        refs.push_back( v[i].at(index) );
    }

    return refs;
}

}


#endif