#ifndef EHDG_SWE_STAB_PARAM_LF_HPP
#define EHDG_SWE_STAB_PARAM_LF_HPP

#include "utilities/sign.hpp"

namespace SWE {
namespace EHDG {
void get_tau_LF(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& tau) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double c  = std::sqrt(Global::g * h);
        double un = u * nx + v * ny;

        tau[gp] = (c + std::abs(un)) * IdentityMatrix<double>(3);
    }
}

void get_dtau_dze_LF(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dtau_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double un      = u * nx + v * ny;
        double dc_dze  = std::sqrt(Global::g / h) / 2.0;
        double dun_dze = -un / h;

        dtau_dze[gp] = (dc_dze + dun_dze * Utilities::sign(un)) * IdentityMatrix<double>(3);
    }
}

void get_dtau_dqx_LF(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dtau_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double un      = u * nx + v * ny;
        double dun_dqx = nx / h;

        dtau_dqx[gp] = dun_dqx * Utilities::sign(un) * IdentityMatrix<double>(3);
    }
}

void get_dtau_dqy_LF(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dtau_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double un      = u * nx + v * ny;
        double dun_dqy = ny / h;

        dtau_dqy[gp] = dun_dqy * Utilities::sign(un) * IdentityMatrix<double>(3);
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
}
}

#endif