#ifndef EHDG_SWE_EDGE_DATA_GLOBAL_HPP
#define EHDG_SWE_EDGE_DATA_GLOBAL_HPP

#include "general_definitions.hpp"
#include <eigen3/Eigen/Dense>

namespace SWE {
namespace EHDG {
struct EdgeGlobal {
    EdgeGlobal() = default;
    EdgeGlobal(const uint ndof, const uint ngp)
        : global_matrix(3 * ndof, 3 * ndof),
          global_rhs(3 * ndof),
          global_delta_q(3 * ndof),
          delta_ze_hat_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_qx_hat_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_qy_hat_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          ze_rhs_kernel_at_gp(ngp),
          qx_rhs_kernel_at_gp(ngp),
          qy_rhs_kernel_at_gp(ngp) {}

    Eigen::MatrixXf global_matrix;
    Eigen::VectorXf global_rhs;
    Eigen::VectorXf global_delta_q;

    std::array<std::vector<double>, 3> delta_ze_hat_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_qx_hat_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_qy_hat_kernel_at_gp;

    std::vector<double> ze_rhs_kernel_at_gp;
    std::vector<double> qx_rhs_kernel_at_gp;
    std::vector<double> qy_rhs_kernel_at_gp;
};
}
}

#endif
