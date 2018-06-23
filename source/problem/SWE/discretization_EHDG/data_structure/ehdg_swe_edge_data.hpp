#ifndef EHDG_SWE_EDGE_DATA_HPP
#define EHDG_SWE_EDGE_DATA_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct EdgeData {
    EdgeData() = default;
    EdgeData(const uint ndof, const uint ngp)
        : ndof(ndof),
          ngp(ngp),
          ze_hat(ndof),
          qx_hat(ndof),
          qy_hat(ndof),
          delta_ze_hat(ndof),
          delta_qx_hat(ndof),
          delta_qy_hat(ndof),
          ze_hat_at_gp(ngp),
          qx_hat_at_gp(ngp),
          qy_hat_at_gp(ngp),
          dT_dze_hat_at_gp(Array3D<double>(3, Array2D<double>(3, std::vector<double>(ngp)))),
          dT_dqx_hat_at_gp(Array3D<double>(3, Array2D<double>(3, std::vector<double>(ngp)))),
          dT_dqy_hat_at_gp(Array3D<double>(3, Array2D<double>(3, std::vector<double>(ngp)))),
          dF_hat_dze_hat_at_gp(Array3D<double>(3, Array2D<double>(2, std::vector<double>(ngp)))),
          dF_hat_dqx_hat_at_gp(Array3D<double>(3, Array2D<double>(2, std::vector<double>(ngp)))),
          dF_hat_dqy_hat_at_gp(Array3D<double>(3, Array2D<double>(2, std::vector<double>(ngp)))) {}

    std::vector<double> ze_hat;
    std::vector<double> qx_hat;
    std::vector<double> qy_hat;

    std::vector<double> delta_ze_hat;
    std::vector<double> delta_qx_hat;
    std::vector<double> delta_qy_hat;

    std::vector<double> ze_hat_at_gp;
    std::vector<double> qx_hat_at_gp;
    std::vector<double> qy_hat_at_gp;

    Array3D<double> dT_dze_hat_at_gp;
    Array3D<double> dT_dqx_hat_at_gp;
    Array3D<double> dT_dqy_hat_at_gp;

    Array3D<double> dF_hat_dze_hat_at_gp;
    Array3D<double> dF_hat_dqx_hat_at_gp;
    Array3D<double> dF_hat_dqy_hat_at_gp;

    uint get_ndof() { return this->ndof; }
    uint get_ngp() { return this->ngp; }

  private:
    uint ndof;
    uint ngp;
};
}
}

#endif
