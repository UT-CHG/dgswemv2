#ifndef IHDG_SWE_DATA_LOCAL_HPP
#define IHDG_SWE_DATA_LOCAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct Local {
    Local() = default;
    Local(const uint ndof, const uint ngp)
        : delta_ze_kernel_at_gp(ngp),
          delta_qx_kernel_at_gp(ngp),
          delta_qy_kernel_at_gp(ngp),
          delta_ze_grad_x_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_ze_grad_y_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_qx_grad_x_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_qx_grad_y_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_qy_grad_x_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          delta_qy_grad_y_kernel_at_gp({std::vector<double>(ngp), std::vector<double>(ngp), std::vector<double>(ngp)}),
          ze_rhs_kernel_at_gp(ngp),
          qx_rhs_kernel_at_gp(ngp),
          qy_rhs_kernel_at_gp(ngp),
          ze_rhs_grad_x_kernel_at_gp(ngp),
          ze_rhs_grad_y_kernel_at_gp(ngp),
          qx_rhs_grad_x_kernel_at_gp(ngp),
          qx_rhs_grad_y_kernel_at_gp(ngp),
          qy_rhs_grad_x_kernel_at_gp(ngp),
          qy_rhs_grad_y_kernel_at_gp(ngp) {}

    std::vector<double> delta_ze_kernel_at_gp;
    std::vector<double> delta_qx_kernel_at_gp;
    std::vector<double> delta_qy_kernel_at_gp;

    std::array<std::vector<double>, 3> delta_ze_grad_x_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_ze_grad_y_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_qx_grad_x_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_qx_grad_y_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_qy_grad_x_kernel_at_gp;
    std::array<std::vector<double>, 3> delta_qy_grad_y_kernel_at_gp;

    std::vector<double> ze_rhs_kernel_at_gp;
    std::vector<double> qx_rhs_kernel_at_gp;
    std::vector<double> qy_rhs_kernel_at_gp;

    std::vector<double> ze_rhs_grad_x_kernel_at_gp;
    std::vector<double> ze_rhs_grad_y_kernel_at_gp;
    std::vector<double> qx_rhs_grad_x_kernel_at_gp;
    std::vector<double> qx_rhs_grad_y_kernel_at_gp;
    std::vector<double> qy_rhs_grad_x_kernel_at_gp;
    std::vector<double> qy_rhs_grad_y_kernel_at_gp;
};
}
}

#endif
