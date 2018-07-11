#ifndef IHDG_SWE_EDGE_DATA_GLOBAL_HPP
#define IHDG_SWE_EDGE_DATA_GLOBAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeGlobal {
    EdgeGlobal() = default;
    EdgeGlobal(const uint ndof_global, const uint ndof_local, const uint ngp)
        : delta_matrix(SWE::n_variables * ndof_global, SWE::n_variables * ndof_local),
          delta_hat_matrix(SWE::n_variables * ndof_global, SWE::n_variables * ndof_global),
          rhs(SWE::n_variables * ndof_global),
          delta_hat_kernel_at_gp(ngp),
          rhs_kernel_at_gp(ngp) {}

    DMatrix<double> delta_matrix;
    DMatrix<double> delta_hat_matrix;
    DVector<double> rhs;

    std::vector<Matrix<double, SWE::n_variables, SWE::n_variables>> delta_kernel_at_gp;
    std::vector<Matrix<double, SWE::n_variables, SWE::n_variables>> delta_hat_kernel_at_gp;
    std::vector<Vector<double, SWE::n_variables>> rhs_kernel_at_gp;
};
}
}

#endif
