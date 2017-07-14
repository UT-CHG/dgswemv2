#ifndef ALMOST_EQUAL_HPP
#define ALMOST_EQUAL_HPP

#include "../general_definitions.hpp"

namespace Utilities {
        constexpr bool almost_equal(double a, double b, double factor = 100)
        {
        //http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
        // the machine epsilon has to be scaled to the magnitude of the values used
        // and multiplied by the desired precision in ULPs (units in the last place)
                return std::abs(a - b) < std::numeric_limits<double>::epsilon() * std::abs(x + y) * factor
                // unless the result is subnormal
                || std::abs(a - b) < std::numeric_limits<double>::min();
        };
}

#endif