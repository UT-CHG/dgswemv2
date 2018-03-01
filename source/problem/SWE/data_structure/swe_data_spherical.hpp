#ifndef SWE_DATA_SPHERICAL_HPP
#define SWE_DATA_SPHERICAL_HPP

#include "../../../general_definitions.hpp"

namespace SWE {
struct Spherical {
    Spherical() = default;
    Spherical(const uint ndof, const uint nbound, const uint ngp_internal, const std::vector<uint>& ngp_boundary) 
      : x(ndof), y(ndof) {
        this->sp_at_gp_internal = std::vector<double>(ngp_internal, 1.0);

        for(uint bound_id=0;bound_id<nbound;bound_id++){
            this->sp_at_gp_boundary.push_back(std::vector<double>(ngp_boundary[bound_id],1.0));
        }
    }

    std::vector<double> x;
    std::vector<double> y;

    std::vector<double> sp_at_gp_internal;
    std::vector<std::vector<double>> sp_at_gp_boundary;    
};
}

#endif