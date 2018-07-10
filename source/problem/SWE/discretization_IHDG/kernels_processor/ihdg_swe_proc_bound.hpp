#ifndef IHDG_SWE_PROC_BOUND_HPP
#define IHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace IHDG {
template <typename BoundaryType>
void Problem::prepare_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage + 1];
    auto& boundary = bound.data.boundary[bound.bound_id];

    bound.ComputeUgp(state.q, boundary.q_at_gp);

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        boundary.aux_at_gp[gp][SWE::Auxiliaries::h] =
            boundary.q_at_gp[gp][SWE::Variables::ze] + boundary.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    /* Compute fluxes at boundary states */
    double nx, ny;
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        nx = bound.surface_normal[gp][GlobalCoord::x];
        ny = bound.surface_normal[gp][GlobalCoord::y];

        u = boundary.q_at_gp[gp][SWE::Variables::qx] / boundary.aux_at_gp[gp][SWE::Auxiliaries::h];
        v = boundary.q_at_gp[gp][SWE::Variables::qy] / boundary.aux_at_gp[gp][SWE::Auxiliaries::h];

        uuh = u * boundary.q_at_gp[gp][SWE::Variables::qx];
        vvh = v * boundary.q_at_gp[gp][SWE::Variables::qy];
        uvh = u * boundary.q_at_gp[gp][SWE::Variables::qy];
        pe  = Global::g * (0.5 * std::pow(boundary.q_at_gp[gp][SWE::Variables::ze], 2) +
                          boundary.q_at_gp[gp][SWE::Variables::ze] * boundary.aux_at_gp[gp][SWE::Auxiliaries::bath]);

        // Fn terms
        boundary.Fn_at_gp[gp][SWE::Variables::ze] =
            boundary.q_at_gp[gp][SWE::Variables::qx] * nx + boundary.q_at_gp[gp][SWE::Variables::qy] * ny;
        boundary.Fn_at_gp[gp][SWE::Variables::qx] = (uuh + pe) * nx + uvh * ny;
        boundary.Fn_at_gp[gp][SWE::Variables::qy] = uvh * nx + (vvh + pe) * ny;

        // dFn/dq terms
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qx) = nx;
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qy) = ny;

        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::ze) =
            (-u * u + Global::g * boundary.aux_at_gp[gp][SWE::Auxiliaries::h]) * nx + -u * v * ny;
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qx) = 2 * u * nx + v * ny;
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qy) = u * ny;

        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::ze) =
            -u * v * nx + (-v * v + Global::g * boundary.aux_at_gp[gp][SWE::Auxiliaries::h]) * ny;
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qx) = v * nx;
        boundary.dF_hat_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qy) = u * nx + 2 * v * ny;
    }
}
}
}

#endif
