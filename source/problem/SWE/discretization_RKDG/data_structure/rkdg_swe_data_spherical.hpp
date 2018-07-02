#ifndef RKDG_SWE_DATA_SPHERICAL_HPP
#define RKDG_SWE_DATA_SPHERICAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct Spherical {
    Spherical() = default;
    Spherical(const uint nnode, const uint nbound, const uint ngp_internal, const std::vector<uint>& ngp_boundary)
        : x_node(nnode), y_node(nnode) {
        this->sp_at_gp_internal = std::vector<double>(ngp_internal, 1.0);

        for (uint bound_id = 0; bound_id < nbound; bound_id++) {
            this->sp_at_gp_boundary.push_back(std::vector<double>(ngp_boundary[bound_id], 1.0));
        }
    }

    std::vector<double> x_node;
    std::vector<double> y_node;

    std::vector<double> sp_at_gp_internal;
    std::vector<std::vector<double>> sp_at_gp_boundary;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void Spherical::serialize(Archive& ar, unsigned) {
    // clang-format off
    ar  & x_node   
        & y_node
        & sp_at_gp_internal
        & sp_at_gp_boundary;
    // clang-format on
}
#endif
}
}

#endif