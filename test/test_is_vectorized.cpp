#include <iostream>
#include <tuple>
#include <utility>

#include "utilities/is_vectorized.hpp"

struct A {
    static constexpr bool is_vectorized = true;
};

struct B {
    static constexpr bool is_vectorized = false;
};

template<typename T>
typename std::enable_if<Utilities::is_vectorized<T>::value>::type function(const T& t) {
    std::cout << "Calling implementation with trait\n";
};

template<typename T>
typename std::enable_if<!Utilities::is_vectorized<T>::value>::type function(const T& t) {
    std::cout << "Calling implementation without trait\n";
};


int main() {

    A a;
    B b;
    using C = std::tuple<>;
    C c;


    std::cout << std::boolalpha << "is_vectorized<A>::value: " << Utilities::is_vectorized<A>::value << '\n'
              << "is_vectorized<B>::value: " << Utilities::is_vectorized<B>::value << '\n'
              << "is_vectorized<C>::value: " << Utilities::is_vectorized<C>::value << '\n';
    function(a);
    function(b);
    function(c);

    return 0;
}