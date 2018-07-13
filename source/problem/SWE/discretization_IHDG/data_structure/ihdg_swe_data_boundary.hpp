#ifndef IHDG_SWE_DATA_BOUNDARY_HPP
#define IHDG_SWE_DATA_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : q_at_gp(ngp),
          aux_at_gp(ngp),
          F_hat_at_gp(ngp),
          dF_hat_dq_at_gp(ngp),
          dF_hat_dq_hat_at_gp(ngp),
          delta_global_kernel_at_gp(ngp) {}

    std::vector<Vector<double, SWE::n_variables>> q_at_gp;
    std::vector<Vector<double, SWE::n_auxiliaries>> aux_at_gp;

    std::vector<Vector<double, SWE::n_variables>> F_hat_at_gp;
    std::vector<Matrix<double, SWE::n_variables, SWE::n_variables>> dF_hat_dq_at_gp;
    std::vector<Matrix<double, SWE::n_variables, SWE::n_variables>> dF_hat_dq_hat_at_gp;
    std::vector<Matrix<double, SWE::n_variables, SWE::n_variables>> delta_global_kernel_at_gp;

    DMatrix<double> delta_global;
    DMatrix<double> delta_hat_local;

    uint ElementID;
    uint EdgeID;
};
}
}

#endif
