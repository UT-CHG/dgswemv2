#ifndef GREENGAUSS_DATA_DERIVATIVE_HPP
#define GREENGAUSS_DATA_DERIVATIVE_HPP

#include "interpolation_data.hpp"

namespace GN {
struct Derivative : GN::Interpolation {
    Derivative() = default;
    Derivative(const uint nvrtx, const uint nbound, const std::vector<uint>& ngp_boundary)
        : GN::Interpolation(nvrtx, nbound, ngp_boundary) {}
};
}

#endif
