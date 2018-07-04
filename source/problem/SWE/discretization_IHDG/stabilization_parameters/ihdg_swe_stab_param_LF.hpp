#ifndef IHDG_SWE_STAB_PARAM_LF_HPP
#define IHDG_SWE_STAB_PARAM_LF_HPP

#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
inline void add_kernel_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

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

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_state.h_hat_at_gp[gp]);

        /* delta kernels */

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

        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] += -2 * tau;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] += -2 * tau;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] += -2 * tau;

        /* RHS kernels */

        edge_global.ze_rhs_kernel_at_gp[gp] += -(tau * (boundary_in.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
                                                 tau * (boundary_ex.ze_at_gp[gp_ex] - edge_state.ze_hat_at_gp[gp]));
        edge_global.qx_rhs_kernel_at_gp[gp] += -(tau * (boundary_in.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
                                                 tau * (boundary_ex.qx_at_gp[gp_ex] - edge_state.qx_hat_at_gp[gp]));
        edge_global.qy_rhs_kernel_at_gp[gp] += -(tau * (boundary_in.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
                                                 tau * (boundary_ex.qy_at_gp[gp_ex] - edge_state.qy_hat_at_gp[gp]));
    }
}

template <typename EdgeDistributedType>
inline void add_kernel_tau_terms_dbound_LF(EdgeDistributedType& edge_dbound) {
    auto& edge_state  = edge_dbound.edge_data.edge_state;
    auto& edge_global = edge_dbound.edge_data.edge_global;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    double ze_ex, qx_ex, qy_ex;
    double ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(
            gp, ze_ex, qx_ex, qy_ex, ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex);

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_dbound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_dbound.boundary.surface_normal[gp][GlobalCoord::y];

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

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_state.h_hat_at_gp[gp]);

        /* delta kernels */

        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] =
            dtau_dze_hat * (boundary.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
            dtau_dze_hat * (ze_ex - edge_state.ze_hat_at_gp[gp]);
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qx][gp] =
            dtau_dze_hat * (boundary.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
            dtau_dze_hat * (qx_ex - edge_state.qx_hat_at_gp[gp]);
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qy][gp] =
            dtau_dze_hat * (boundary.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
            dtau_dze_hat * (qy_ex - edge_state.qy_hat_at_gp[gp]);

        edge_global.delta_qx_hat_kernel_at_gp[Variables::ze][gp] =
            dtau_dqx_hat * (boundary.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
            dtau_dqx_hat * (ze_ex - edge_state.ze_hat_at_gp[gp]);
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] =
            dtau_dqx_hat * (boundary.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
            dtau_dqx_hat * (qx_ex - edge_state.qx_hat_at_gp[gp]);
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qy][gp] =
            dtau_dqx_hat * (boundary.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
            dtau_dqx_hat * (qy_ex - edge_state.qy_hat_at_gp[gp]);

        edge_global.delta_qy_hat_kernel_at_gp[Variables::ze][gp] =
            dtau_dqy_hat * (boundary.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
            dtau_dqy_hat * (ze_ex - edge_state.ze_hat_at_gp[gp]);
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qx][gp] =
            dtau_dqy_hat * (boundary.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
            dtau_dqy_hat * (qx_ex - edge_state.qx_hat_at_gp[gp]);
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] =
            dtau_dqy_hat * (boundary.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
            dtau_dqy_hat * (qy_ex - edge_state.qy_hat_at_gp[gp]);

        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] += -2 * tau;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] += -2 * tau;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] += -2 * tau;

        /* RHS kernels */

        edge_global.ze_rhs_kernel_at_gp[gp] += -(tau * (boundary.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]) +
                                                 tau * (ze_ex - edge_state.ze_hat_at_gp[gp]));
        edge_global.qx_rhs_kernel_at_gp[gp] += -(tau * (boundary.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]) +
                                                 tau * (qx_ex - edge_state.qx_hat_at_gp[gp]));
        edge_global.qy_rhs_kernel_at_gp[gp] += -(tau * (boundary.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]) +
                                                 tau * (qy_ex - edge_state.qy_hat_at_gp[gp]));
    }
}

template <typename EdgeInterfaceType>
inline void add_flux_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_state = edge_int.edge_data.edge_state;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

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

template <typename EdgeBoundaryType>
inline void add_flux_tau_terms_bound_LF(EdgeBoundaryType& edge_bound) {
    auto& edge_state = edge_bound.edge_data.edge_state;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_state.h_hat_at_gp[gp]);

        boundary.ze_numerical_flux_at_gp[gp] += tau * (boundary.ze_at_gp[gp] - edge_state.ze_hat_at_gp[gp]);
        boundary.qx_numerical_flux_at_gp[gp] += tau * (boundary.qx_at_gp[gp] - edge_state.qx_hat_at_gp[gp]);
        boundary.qy_numerical_flux_at_gp[gp] += tau * (boundary.qy_at_gp[gp] - edge_state.qy_hat_at_gp[gp]);
    }
}
}
}

#endif