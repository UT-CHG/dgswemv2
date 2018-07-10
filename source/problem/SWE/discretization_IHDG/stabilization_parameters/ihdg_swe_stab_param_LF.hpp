#ifndef IHDG_SWE_STAB_PARAM_LF_HPP
#define IHDG_SWE_STAB_PARAM_LF_HPP

#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
inline void add_dF_hat_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

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

        dtau_dze_hat = std::sqrt(Global::g / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]) / 2.0 -
                       sgn * un_hat / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dqx_hat = sgn * nx / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dqy_hat = sgn * ny / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        Vector<double, SWE::n_variables> del_u_in = boundary_in.q_at_gp[gp] - edge_internal.q_hat_at_gp[gp];
        Vector<double, SWE::n_variables> del_u_ex = boundary_ex.q_at_gp[gp_ex] - edge_internal.q_hat_at_gp[gp];

        // dF_hat_dq tau terms
        boundary_in.dF_hat_dq_at_gp[gp]    += tau * blaze::IdentityMatrix<double>(SWE::n_variables);
        boundary_ex.dF_hat_dq_at_gp[gp_ex] += tau * blaze::IdentityMatrix<double>(SWE::n_variables);

        // dF_hat_dq_hat tau terms
        blaze::column<SWE::Variables::ze>(boundary_in.dF_hat_dq_hat_at_gp[gp]) = dtau_dze_hat * del_u_in;
        blaze::column<SWE::Variables::qx>(boundary_in.dF_hat_dq_hat_at_gp[gp]) = dtau_dqx_hat * del_u_in;
        blaze::column<SWE::Variables::qy>(boundary_in.dF_hat_dq_hat_at_gp[gp]) = dtau_dqy_hat * del_u_in;

        boundary_in.dF_hat_dq_hat_at_gp[gp] += -tau * blaze::IdentityMatrix<double>(SWE::n_variables);

        blaze::column<SWE::Variables::ze>(boundary_ex.dF_hat_dq_hat_at_gp[gp_ex]) = dtau_dze_hat * del_u_ex;
        blaze::column<SWE::Variables::qx>(boundary_ex.dF_hat_dq_hat_at_gp[gp_ex]) = dtau_dqx_hat * del_u_ex;
        blaze::column<SWE::Variables::qy>(boundary_ex.dF_hat_dq_hat_at_gp[gp_ex]) = dtau_dqy_hat * del_u_ex;

        boundary_ex.dF_hat_dq_hat_at_gp[gp] += -tau * blaze::IdentityMatrix<double>(SWE::n_variables);
    }
}

template <typename EdgeBoundaryType>
inline void add_dF_hat_tau_terms_boundary_LF(EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        if (0.0 < un_hat) {
            sgn = 1.0;
        } else if (un_hat < 0.0) {
            sgn = -1.0;
        } else {
            sgn = 0.0;
        }

        dtau_dze_hat = std::sqrt(Global::g / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]) / 2.0 -
                       sgn * un_hat / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dqx_hat = sgn * nx / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dqy_hat = sgn * ny / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        Vector<double, SWE::n_variables> del_u = boundary.q_at_gp[gp] - edge_internal.q_hat_at_gp[gp];

        // dF_hat_dq tau terms
        boundary.dF_hat_dq_at_gp[gp]    += tau * blaze::IdentityMatrix<double>(SWE::n_variables);

        // dF_hat_dq_hat tau terms
        blaze::column<SWE::Variables::ze>(boundary.dF_hat_dq_hat_at_gp[gp]) = dtau_dze_hat * del_u;
        blaze::column<SWE::Variables::qx>(boundary.dF_hat_dq_hat_at_gp[gp]) = dtau_dqx_hat * del_u;
        blaze::column<SWE::Variables::qy>(boundary.dF_hat_dq_hat_at_gp[gp]) = dtau_dqy_hat * del_u;

        boundary.dF_hat_dq_hat_at_gp[gp] += -tau * blaze::IdentityMatrix<double>(SWE::n_variables);
    }
}

template <typename EdgeDistributedType>
inline void add_kernel_tau_terms_dbound_LF(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;
    auto& edge_global   = edge_dbound.edge_data.edge_global;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    Vector<double, SWE::n_variables> q_ex;
    Vector<double, SWE::n_variables> Fn_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(gp, q_ex, Fn_ex);

        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

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

        dtau_dze_hat = std::sqrt(Global::g / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]) / 2.0 -
                       sgn * un_hat / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dqx_hat = sgn * nx / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dqy_hat = sgn * ny / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        Vector<double, SWE::n_variables> del_u = boundary.q_at_gp[gp] + q_ex - 2.0 * edge_internal.q_hat_at_gp[gp];

        /* delta kernels */

        blaze::column<SWE::Variables::ze>(edge_global.delta_hat_kernel_at_gp[gp]) = dtau_dze_hat * del_u;
        blaze::column<SWE::Variables::qx>(edge_global.delta_hat_kernel_at_gp[gp]) = dtau_dqx_hat * del_u;
        blaze::column<SWE::Variables::qy>(edge_global.delta_hat_kernel_at_gp[gp]) = dtau_dqy_hat * del_u;

        edge_global.delta_hat_kernel_at_gp[gp] += -2 * tau * blaze::IdentityMatrix<double>(SWE::n_variables);

        /* RHS kernels */

        edge_global.rhs_kernel_at_gp[gp] += -tau * del_u;
    }
}

template <typename EdgeInterfaceType>
inline void add_F_hat_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        boundary_in.F_hat_at_gp[gp] += tau * (boundary_in.q_at_gp[gp] - edge_internal.q_hat_at_gp[gp]);
        boundary_ex.F_hat_at_gp[gp_ex] += tau * (boundary_ex.q_at_gp[gp_ex] - edge_internal.q_hat_at_gp[gp]);
    }
}

template <typename EdgeBoundaryType>
inline void add_F_hat_tau_terms_bound_LF(EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        boundary.F_hat_at_gp[gp] += tau * (boundary.q_at_gp[gp] - edge_internal.q_hat_at_gp[gp]);
    }
}
}
}

#endif