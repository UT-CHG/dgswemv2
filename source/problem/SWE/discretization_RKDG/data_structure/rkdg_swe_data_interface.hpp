#ifndef RKDG_SWE_DATA_INTERFACE_HPP
#define RKDG_SWE_DATA_INTERFACE_HPP

namespace SWE {
namespace RKDG {

struct InterfaceData {
    InterfaceData()=default;
    InterfaceData(const uint ninterfaces, const uint ngp) {
        q_in_at_gp.fill(DynMatrix<double>(ninterfaces, ngp));
        q_ex_at_gp.fill(DynMatrix<double>(ninterfaces, ngp));

        aux_at_gp.fill(DynMatrix<double>(ninterfaces, ngp));

        F_hat_at_gp.fill(DynMatrix<double>(ninterfaces, ngp));
    }

    std::array<DynMatrix<double>, SWE::n_variables> q_in_at_gp;
    std::array<DynMatrix<double>, SWE::n_variables> q_ex_at_gp;

    std::array<DynMatrix<double>, SWE::n_auxiliaries> aux_at_gp;

    std::array<DynMatrix<double>, SWE::n_variables> F_hat_at_gp;
};
}
}

#endif