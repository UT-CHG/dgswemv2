#ifndef SWE_INITIAL_CONDITION_FUNCTIONS_HPP
#define SWE_INITIAL_CONDITION_FUNCTIONS_HPP

#include "../../../general_definitions.hpp"

namespace SWE {
inline double ic_ze(const double t, const Point<2>& pt) { return 0; }

inline double ic_qx(const double t, const Point<2>& pt) { return 0; }

inline double ic_qy(const double t, const Point<2>& pt) { return 0; }
}

#endif
