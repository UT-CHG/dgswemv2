#ifndef EHDG_SWE_EDGE_DATA_INTERNAL_HPP
#define EHDG_SWE_EDGE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp) : q_avg_at_gp(ngp), q_hat_at_gp(ngp), h_hat_at_gp(ngp) {}

    std::vector<Vector<double, SWE::n_variables>> q_avg_at_gp;
    std::vector<Vector<double, SWE::n_variables>> q_hat_at_gp;
    std::vector<double> h_hat_at_gp;
};
}
}

#endif
