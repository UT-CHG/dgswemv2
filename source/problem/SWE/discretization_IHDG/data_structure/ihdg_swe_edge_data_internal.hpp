#ifndef IHDG_SWE_EDGE_DATA_INTERNAL_HPP
#define IHDG_SWE_EDGE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : q_hat_at_gp(SWE::n_variables, ngp),
          aux_hat_at_gp(SWE::n_auxiliaries, ngp),
          q_init_at_gp(SWE::n_variables, ngp),
          delta_hat_global_kernel_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          rhs_global_kernel_at_gp(SWE::n_variables, ngp) {}

    DynMatrix<double> q_hat_at_gp;
    DynMatrix<double> aux_hat_at_gp;
    DynMatrix<double> q_init_at_gp;

    DynMatrix<double> delta_hat_global_kernel_at_gp;
    DynMatrix<double> rhs_global_kernel_at_gp;

    DynMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;
};
}
}

#endif
