#ifndef EHDG_SWE_DATA_BOUNDARY_HPP
#define EHDG_SWE_DATA_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp) : q_at_gp(ngp), aux_at_gp(ngp), Fn_at_gp(ngp), F_hat_at_gp(ngp) {}

    std::vector<Vector<double, SWE::n_variables>> q_at_gp;
    std::vector<Vector<double, SWE::n_auxiliaries>> aux_at_gp;

    std::vector<Vector<double, SWE::n_variables>> Fn_at_gp;
    std::vector<Vector<double, SWE::n_variables>> F_hat_at_gp;
};
}
}

#endif
