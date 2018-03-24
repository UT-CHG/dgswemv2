#include "utilities/tuple_helpers.hpp"

#include <iostream>
#include <vector>

int main() {

    std::vector<double> foo{0., 1.1, 5.0};
    std::string         _str = "bar";

    auto tup = std::make_tuple(5, _str, foo);

    std::cout << "Index of foo is " << Utilities::index<decltype(foo), decltype(tup)>::value << std::endl;

    static_assert(Utilities::has_type<int, decltype(tup)>::value, "Error in has_type");
    static_assert(Utilities::has_type<std::string, decltype(tup)>::value, "Error in has_type");
    static_assert(Utilities::has_type<int, decltype(tup)>::value, "Error in has_type");
};
