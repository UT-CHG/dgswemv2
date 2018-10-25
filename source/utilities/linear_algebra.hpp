#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#define USE_EIGEN
//#define USE_BLAZE

#ifdef USE_BLAZE
#include "utilities/linear_algebra/use_blaze.hpp"
#endif

#ifdef USE_EIGEN
#include "utilities/linear_algebra/use_eigen.hpp"
#endif


//The following are STL containers with aligned allocators.
//These should be used whenever the template parameter is
//a static or Hybrid vector type or contains is a class
//which contains Hybrid or static vector types
template<typename T>
using AlignedVector = std::vector<T,AlignedAllocator<T>>;

template<typename Key, typename T, typename Compare = std::less<Key> >
using AlignedMap = std::map<Key,T,Compare,AlignedAllocator<std::pair<Key,T>>>;
#endif