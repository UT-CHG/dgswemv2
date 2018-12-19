#include <iostream>
#include <tuple>
#include <utility>

#include "utilities/is_vectorized.hpp"

struct A {
    template <typename Arg>
    static constexpr bool is_vectorized() { return true; }
};

struct B {
    template <typename Arg>
    static constexpr bool is_vectorized() { return false; }
};

struct C {

    template <typename Arg>
    static constexpr bool is_vectorized() { return false; }
};

template <>
constexpr bool C::is_vectorized<double>() { return true; }

template<typename T, typename AccessorType = std::tuple<> >
typename std::enable_if<Utilities::is_vectorized<T, AccessorType>::value>::type function(const T& t) {
    std::cout << "Calling implementation with trait\n";
};

template<typename T, typename AccessorType = std::tuple<> >
typename std::enable_if<!Utilities::is_vectorized<T, AccessorType>::value>::type function(const T& t) {
    std::cout << "Calling implementation without trait\n";
};


int main() {

    A a;
    B b;
    C c;
    using D = std::tuple<>;
    D d;
    auto e = [](auto& var) {};

    using default_t = std::tuple<>;

    std::cout << std::boolalpha
              << "is_vectorized<A>::value: " << Utilities::is_vectorized<A>::value << '\n'
              << "is_vectorized<B>::value: " << Utilities::is_vectorized<B>::value << '\n'
              << "is_vectorized<C>::value: " << Utilities::is_vectorized<C>::value << '\n'
              << "is_vectorized<C,double>::value " << Utilities::is_vectorized<C,double>::value << '\n'
              << "is_vectorized<D>::value " << Utilities::is_vectorized<D>::value << '\n'
              << "is_vectorized<decltype(d)>::value: " << Utilities::is_vectorized<decltype(e)>::value << '\n';

    function(a);
    function(b);
    function(c);
    function<C,double>(c);
    function(d);
    function(e);

    return ( Utilities::is_vectorized<A>::value != true ) ||
        (Utilities::is_vectorized<B>::value != false ) ||
        (Utilities::is_vectorized<C>::value != false) ||
        (Utilities::is_vectorized<C, double>::value != true) ||
        (Utilities::is_vectorized<D>::value != false) ||
        (Utilities::is_vectorized<decltype(e)>::value != false );
}