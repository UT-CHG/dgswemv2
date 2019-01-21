#ifndef RKDG_SWE_DATA_INTERFACE_HPP
#define RKDG_SWE_DATA_INTERFACE_HPP

namespace SWE {
namespace RKDG {

struct InterfaceData {
    InterfaceData()=default;
    InterfaceData(const uint ninterfaces, const uint ngp) {
      q_in_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(ninterfaces, ngp));
        q_ex_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(ninterfaces, ngp));

        aux_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(ninterfaces, ngp));

        F_hat_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(ninterfaces, ngp));
    }

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> q_in_at_gp;
    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> q_ex_at_gp;

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp;

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> F_hat_at_gp;
};
}
}

#endif