#ifndef GREENGAUSS_DATA_DERIVATIVE_HPP
#define GREENGAUSS_DATA_DERIVATIVE_HPP

#include "reconstruction_data.hpp"

namespace GN {
struct Derivative : GN::Reconstruction {
    Derivative() = default;
    Derivative(const uint nvrtx, const uint nbound, const std::vector<uint>& ngp_boundary)
        : GN::Reconstruction(nvrtx, nbound, ngp_boundary) {}
};
}

#endif
