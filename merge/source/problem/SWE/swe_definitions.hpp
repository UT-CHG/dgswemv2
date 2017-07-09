#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

namespace SWE {
	namespace Global {
		static constexpr double g = 9.81;
		static constexpr double Cf = 0.0025;
	}

	enum BoundaryConditions : unsigned char {
		land = 0,
		tidal = 1,
		internal = 255
	};
}

#endif
