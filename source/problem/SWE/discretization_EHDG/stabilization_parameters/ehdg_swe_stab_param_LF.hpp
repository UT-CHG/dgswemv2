#ifndef EHDG_SWE_STAB_PARAM_LF_HPP
#define EHDG_SWE_STAB_PARAM_LF_HPP

#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
namespace EHDG {
template <typename EdgeInternalType>
inline void add_delta_kernel_tau_terms_LF(EdgeInternalType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.data_in.boundary[edge_int.bound_id_in];
    auto& boundary_ex = edge_int.data_ex.boundary[edge_int.bound_id_ex];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_int.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.surface_normal_in[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        if (0.0 < un_hat) {
            sgn = 1.0;
        } else if (un_hat < 0.0) {
            sgn = -1.0;
        } else {
            sgn = 0.0;
        }

        dtau_dze_hat =
            std::sqrt(Global::g / edge_state.h_hat_at_gp[gp]) / 2.0 - sgn * un_hat / edge_state.h_hat_at_gp[gp];
        dtau_dqx_hat = sgn * nx / edge_state.h_hat_at_gp[gp];
        dtau_dqy_hat = sgn * ny / edge_state.h_hat_at_gp[gp];

        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] =
            dtau_dze_hat * (boundary_in.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
            dtau_dze_hat * (boundary_ex.ze_at_gp[gp_ex] - edge_state.ze_hat_at_gp[gp]);
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qx][gp] =
            dtau_dze_hat * (boundary_in.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
            dtau_dze_hat * (boundary_ex.qx_at_gp[gp_ex] - edge_state.qx_hat_at_gp[gp]);
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qy][gp] =
            dtau_dze_hat * (boundary_in.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
            dtau_dze_hat * (boundary_ex.qy_at_gp[gp_ex] - edge_state.qy_hat_at_gp[gp]);

        edge_global.delta_qx_hat_kernel_at_gp[Variables::ze][gp] =
            dtau_dqx_hat * (boundary_in.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
            dtau_dqx_hat * (boundary_ex.ze_at_gp[gp_ex] - edge_state.ze_hat_at_gp[gp]);
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] =
            dtau_dqx_hat * (boundary_in.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
            dtau_dqx_hat * (boundary_ex.qx_at_gp[gp_ex] - edge_state.qx_hat_at_gp[gp]);
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qy][gp] =
            dtau_dqx_hat * (boundary_in.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
            dtau_dqx_hat * (boundary_ex.qy_at_gp[gp_ex] - edge_state.qy_hat_at_gp[gp]);

        edge_global.delta_qy_hat_kernel_at_gp[Variables::ze][gp] =
            dtau_dqy_hat * (boundary_in.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
            dtau_dqy_hat * (boundary_ex.ze_at_gp[gp_ex] - edge_state.ze_hat_at_gp[gp]);
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qx][gp] =
            dtau_dqy_hat * (boundary_in.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
            dtau_dqy_hat * (boundary_ex.qx_at_gp[gp_ex] - edge_state.qx_hat_at_gp[gp]);
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] =
            dtau_dqy_hat * (boundary_in.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
            dtau_dqy_hat * (boundary_ex.qy_at_gp[gp_ex] - edge_state.qy_hat_at_gp[gp]);

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_state.h_hat_at_gp[gp]);

        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] += -2 * tau;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] += -2 * tau;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] += -2 * tau;
    }
}

template <typename EdgeInternalType>
inline void add_rhs_kernel_tau_terms_LF(EdgeInternalType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.data_in.boundary[edge_int.bound_id_in];
    auto& boundary_ex = edge_int.data_ex.boundary[edge_int.bound_id_ex];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_int.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.surface_normal_in[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_state.h_hat_at_gp[gp]);

        edge_global.ze_rhs_kernel_at_gp[gp] += -(tau * (boundary_in.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
                                                 tau * (boundary_ex.ze_at_gp[gp_ex] - edge_state.ze_hat_at_gp[gp]));
        edge_global.qx_rhs_kernel_at_gp[gp] += -(tau * (boundary_in.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
                                                 tau * (boundary_ex.qx_at_gp[gp_ex] - edge_state.qx_hat_at_gp[gp]));
        edge_global.qy_rhs_kernel_at_gp[gp] += -(tau * (boundary_in.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
                                                 tau * (boundary_ex.qy_at_gp[gp_ex] - edge_state.qy_hat_at_gp[gp]));
    }
}

template <typename EdgeInternalType>
inline void add_flux_tau_terms_LF(EdgeInternalType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.data_in.boundary[edge_int.bound_id_in];
    auto& boundary_ex = edge_int.data_ex.boundary[edge_int.bound_id_ex];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_int.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.surface_normal_in[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_state.h_hat_at_gp[gp]);

        boundary_in.ze_numerical_flux_at_gp[gp] += tau * (boundary_in.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]);
        boundary_in.qx_numerical_flux_at_gp[gp] += tau * (boundary_in.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]);
        boundary_in.qy_numerical_flux_at_gp[gp] += tau * (boundary_in.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]);

        boundary_ex.ze_numerical_flux_at_gp[gp_ex] += tau * (boundary_ex.ze_at_gp[gp_ex] - edge_state.ze_hat_at_gp[gp]);
        boundary_ex.qx_numerical_flux_at_gp[gp_ex] += tau * (boundary_ex.qx_at_gp[gp_ex] - edge_state.qx_hat_at_gp[gp]);
        boundary_ex.qy_numerical_flux_at_gp[gp_ex] += tau * (boundary_ex.qy_at_gp[gp_ex] - edge_state.qy_hat_at_gp[gp]);
    }
}
}
}

#endif