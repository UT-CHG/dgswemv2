#include "utilities/is_defined.hpp"

#include <iostream>

struct NotDefined;

struct Defined {};

int main() {
    std::cout << "is_defined<NotDefined>: " << std::boolalpha << Utilities::is_defined<NotDefined>::value
              << " (should be false)\n"
              << "is_defined<Defined>: " << std::boolalpha << Utilities::is_defined<Defined>::value
              << " (should be true)" << std::endl;

    return !(!Utilities::is_defined<NotDefined>::value && Utilities::is_defined<Defined>::value);
}