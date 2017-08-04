#ifndef ALMOST_EQUAL_HPP
#define ALMOST_EQUAL_HPP

#include "../general_definitions.hpp"

namespace Utilities {
	constexpr bool almost_equal(double a, double b, double factor = 100) {
		if (std::max(std::abs(a), std::abs(b)) < std::numeric_limits<double>::epsilon() * factor) {
			return true;
		}

		return std::abs(a - b) < (std::max(std::abs(a), std::abs(b)) *
			std::numeric_limits<double>::epsilon() * factor);
	};
}

#endif