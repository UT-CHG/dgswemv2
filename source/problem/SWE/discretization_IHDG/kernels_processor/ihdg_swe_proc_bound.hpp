#ifndef IHDG_SWE_PROC_BOUND_HPP
#define IHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace IHDG {
template <typename BoundaryType>
void Problem::local_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage + 1];
    auto& boundary = bound.data.boundary[bound.bound_id];

    bound.ComputeUgp(state.q, boundary.q_at_gp);

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        boundary.aux_at_gp(SWE::Auxiliaries::h, gp) =
            boundary.q_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
    }

    /* Compute fluxes at boundary states */
    double nx, ny;
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        nx = bound.surface_normal[gp][GlobalCoord::x];
        ny = bound.surface_normal[gp][GlobalCoord::y];

        u = boundary.q_at_gp(SWE::Variables::qx, gp) / boundary.aux_at_gp(SWE::Auxiliaries::h, gp);
        v = boundary.q_at_gp(SWE::Variables::qy, gp) / boundary.aux_at_gp(SWE::Auxiliaries::h, gp);

        uuh = u * boundary.q_at_gp(SWE::Variables::qx, gp);
        vvh = v * boundary.q_at_gp(SWE::Variables::qy, gp);
        uvh = u * boundary.q_at_gp(SWE::Variables::qy, gp);
        pe  = Global::g * (0.5 * std::pow(boundary.q_at_gp(SWE::Variables::ze, gp), 2) +
                          boundary.q_at_gp(SWE::Variables::ze, gp) * boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));

        // Fn terms
        boundary.F_hat_at_gp(SWE::Variables::ze, gp) =
            boundary.q_at_gp(SWE::Variables::qx, gp) * nx + boundary.q_at_gp(SWE::Variables::qy, gp) * ny;
        boundary.F_hat_at_gp(SWE::Variables::qx, gp) = (uuh + pe) * nx + uvh * ny;
        boundary.F_hat_at_gp(SWE::Variables::qy, gp) = uvh * nx + (vvh + pe) * ny;

        // dFn/dq terms
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::ze_ze] = 0.0;
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::ze_qx] = nx;
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::ze_qy] = ny;

        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::qx_ze] =
            (-u * u + Global::g * boundary.aux_at_gp(SWE::Auxiliaries::h, gp)) * nx - u * v * ny;
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::qx_qx] = 2 * u * nx + v * ny;
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::qx_qy] = u * ny;

        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::qy_ze] =
            -u * v * nx + (-v * v + Global::g * boundary.aux_at_gp(SWE::Auxiliaries::h, gp)) * ny;
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::qy_qx] = v * nx;
        boundary.dF_hat_dq_at_gp[gp][JacobianVariables::qy_qy] = u * nx + 2 * v * ny;
    }
}
}
}

#endif
