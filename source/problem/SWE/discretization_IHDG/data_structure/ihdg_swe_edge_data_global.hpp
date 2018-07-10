#ifndef IHDG_SWE_EDGE_DATA_GLOBAL_HPP
#define IHDG_SWE_EDGE_DATA_GLOBAL_HPP

#include "general_definitions.hpp"
#include <eigen3/Eigen/Dense>

namespace SWE {
namespace IHDG {
struct EdgeGlobal {
    EdgeGlobal() = default;
    EdgeGlobal(const uint ndof, const uint ngp)
        : global_matrix(3 * ndof, 3 * ndof),
          global_rhs(3 * ndof),
          global_delta_q(3 * ndof),
          delta_hat_kernel_at_gp(ngp),
          rhs_kernel_at_gp(ngp) {}

    Eigen::MatrixXf global_matrix;
    Eigen::VectorXf global_rhs;
    Eigen::VectorXf global_delta_q;

    std::vector<Matrix<double, SWE::n_variables, SWE::n_variables>> delta_hat_kernel_at_gp;
    std::vector<Vector<double, SWE::n_variables>> rhs_kernel_at_gp;
};
}
}

#endif
