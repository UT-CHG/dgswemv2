#ifndef ALMOST_EQUAL_HPP
#define ALMOST_EQUAL_HPP

#include "general_definitions.hpp"

namespace Utilities {
/**
 * Determine whether two floating point numbers are approximately equal.
 * Determine whether two numbers are within factor of machine precision for double.
 *
 * @param a
 * @param b
 * @param factor Factor by which machine precision will be multiplied for determining equality.
 * @return whether or not the numbers are approximately equal.
 */
constexpr bool almost_equal(double a, double b, double factor = 100) {
    // implementation from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
    return std::abs(a - b) <= std::numeric_limits<double>::epsilon() * std::abs(a + b) * factor ||
           std::abs(a - b) < std::numeric_limits<double>::min();
};
}

#endif