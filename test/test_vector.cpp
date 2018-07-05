#include "utilities/vector.hpp"
#include "general_definitions.hpp"

int main() {
    Utilities::LinearAlgebra::Vector<double> vec(3);

    vec[0] = 3.5;
    vec[1] = 4.4;
    vec[2] = 5.3;

    std::cout << vec[0] << vec[1] << vec[2] << std::endl;

    Utilities::LinearAlgebra::Vector<double> vec2(3);

    vec2 = 3.0 * vec + vec - 2.0 * vec;

    std::cout << vec2[0] << vec2[1] << vec2[2] << std::endl;
};
