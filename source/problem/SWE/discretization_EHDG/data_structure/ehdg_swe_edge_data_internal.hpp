#ifndef EHDG_SWE_EDGE_DATA_INTERNAL_HPP
#define EHDG_SWE_EDGE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : q_hat_at_gp(ngp),
          aux_hat_at_gp(ngp),
          q_init_at_gp(ngp),
          delta_hat_global_kernel_at_gp(ngp),
          rhs_global_kernel_at_gp(ngp) {}

    std::vector<StatVector<double, SWE::n_variables>> q_hat_at_gp;
    std::vector<StatVector<double, SWE::n_auxiliaries>> aux_hat_at_gp;

    std::vector<StatVector<double, SWE::n_variables>> q_init_at_gp;

    std::vector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> delta_hat_global_kernel_at_gp;
    std::vector<StatVector<double, SWE::n_variables>> rhs_global_kernel_at_gp;

    DynMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;
};
}
}

#endif
