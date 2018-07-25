#ifndef EHDG_SWE_STAB_PARAM_LF_HPP
#define EHDG_SWE_STAB_PARAM_LF_HPP

#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
namespace EHDG {
template <typename EdgeInterfaceType>
inline void add_kernel_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau;
    double sgn;

    StatVector<double, SWE::n_variables> del_q;
    StatVector<double, SWE::n_variables> dtau_dq_hat;

    StatVector<double, SWE::n_variables * SWE::n_variables> dtau_delq;
    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qx, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        v_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qy, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        nx = edge_int.interface.surface_normal_in(GlobalCoord::x, gp);
        ny = edge_int.interface.surface_normal_in(GlobalCoord::y, gp);

        un_hat = u_hat * nx + v_hat * ny;

        if (0.0 < un_hat) {
            sgn = 1.0;
        } else if (un_hat < 0.0) {
            sgn = -1.0;
        } else {
            sgn = 0.0;
        }

        dtau_dq_hat[SWE::Variables::ze] =
            std::sqrt(Global::g / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp)) / 2.0 -
            sgn * un_hat / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        dtau_dq_hat[SWE::Variables::qx] = sgn * nx / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        dtau_dq_hat[SWE::Variables::qy] = sgn * ny / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp));

        del_q = column(boundary_in.q_at_gp, gp) + column(boundary_ex.q_at_gp, gp_ex) -
                2.0 * column(edge_internal.q_hat_at_gp, gp);

        /* delta kernels */
        subvector(dtau_delq, SWE::JacobianVariables::ze_ze, SWE::n_variables) = dtau_dq_hat * del_q[SWE::Variables::ze];
        subvector(dtau_delq, SWE::JacobianVariables::qx_ze, SWE::n_variables) = dtau_dq_hat * del_q[SWE::Variables::qx];
        subvector(dtau_delq, SWE::JacobianVariables::qy_ze, SWE::n_variables) = dtau_dq_hat * del_q[SWE::Variables::qy];

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = dtau_delq - 2.0 * tau * I_vector;
        /* RHS kernels */

        column(edge_internal.rhs_global_kernel_at_gp, gp) += tau * del_q;
    }
}

template <typename EdgeDistributedType>
inline void add_kernel_tau_terms_dbound_LF(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau;
    double sgn;

    StatVector<double, SWE::n_variables> del_q;
    StatVector<double, SWE::n_variables> dtau_dq_hat;

    StatVector<double, SWE::n_variables * SWE::n_variables> dtau_delq;
    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        u_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qx, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        v_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qy, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        nx = edge_dbound.boundary.surface_normal(GlobalCoord::x, gp);
        ny = edge_dbound.boundary.surface_normal(GlobalCoord::y, gp);

        un_hat = u_hat * nx + v_hat * ny;

        if (0.0 < un_hat) {
            sgn = 1.0;
        } else if (un_hat < 0.0) {
            sgn = -1.0;
        } else {
            sgn = 0.0;
        }

        dtau_dq_hat[SWE::Variables::ze] =
            std::sqrt(Global::g / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp)) / 2.0 -
            sgn * un_hat / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        dtau_dq_hat[SWE::Variables::qx] = sgn * nx / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        dtau_dq_hat[SWE::Variables::qy] = sgn * ny / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp));

        del_q = column(boundary.q_at_gp, gp) + column(edge_dbound.boundary.boundary_condition.q_ex, gp) -
                2.0 * column(edge_internal.q_hat_at_gp, gp);

        /* delta kernels */

        subvector(dtau_delq, SWE::JacobianVariables::ze_ze, SWE::n_variables) = dtau_dq_hat * del_q[SWE::Variables::ze];
        subvector(dtau_delq, SWE::JacobianVariables::qx_ze, SWE::n_variables) = dtau_dq_hat * del_q[SWE::Variables::qx];
        subvector(dtau_delq, SWE::JacobianVariables::qy_ze, SWE::n_variables) = dtau_dq_hat * del_q[SWE::Variables::qy];

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = dtau_delq - 2.0 * tau * I_vector;

        /* RHS kernels */

        column(edge_internal.rhs_global_kernel_at_gp, gp) += tau * del_q;
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

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qx, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        v_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qy, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        nx = edge_int.interface.surface_normal_in(GlobalCoord::x, gp);
        ny = edge_int.interface.surface_normal_in(GlobalCoord::y, gp);

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp));

        column(boundary_in.F_hat_at_gp, gp) +=
            tau * (column(boundary_in.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
        column(boundary_ex.F_hat_at_gp, gp_ex) +=
            tau * (column(boundary_ex.q_at_gp, gp_ex) - column(edge_internal.q_hat_at_gp, gp));
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
            edge_internal.q_hat_at_gp(SWE::Variables::qx, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        v_hat =
            edge_internal.q_hat_at_gp(SWE::Variables::qy, gp) / edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        nx = edge_bound.boundary.surface_normal(GlobalCoord::x, gp);
        ny = edge_bound.boundary.surface_normal(GlobalCoord::y, gp);

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp));

        column(boundary.F_hat_at_gp, gp) +=
            tau * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
    }
}
}
}

#endif