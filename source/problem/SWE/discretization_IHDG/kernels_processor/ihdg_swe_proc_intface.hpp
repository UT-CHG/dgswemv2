#ifndef IHDG_SWE_PROC_INTFACE_HPP
#define IHDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace IHDG {
template <typename InterfaceType>
void Problem::local_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage + 1];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage + 1];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    intface.ComputeUgpIN(state_in.q, boundary_in.q_at_gp);
    intface.ComputeUgpEX(state_ex.q, boundary_ex.q_at_gp);

    uint gp_ex;
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        gp_ex = intface.data_in.get_ngp_boundary(intface.bound_id_in) - gp - 1;

        boundary_in.aux_at_gp[gp][SWE::Auxiliaries::h] =
            boundary_in.q_at_gp[gp][SWE::Variables::ze] + boundary_in.aux_at_gp[gp][SWE::Auxiliaries::bath];
        boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::h] =
            boundary_ex.q_at_gp[gp_ex][SWE::Variables::ze] + boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::bath];
    }

    /* Compute fluxes at boundary states */
    double nx_in, ny_in;
    double u_in, v_in;
    double uuh_in, vvh_in, uvh_in, pe_in;

    double nx_ex, ny_ex;
    double u_ex, v_ex;
    double uuh_ex, vvh_ex, uvh_ex, pe_ex;
  
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        gp_ex = intface.data_in.get_ngp_boundary(intface.bound_id_in) - gp - 1;

        /* IN State */

        nx_in = intface.surface_normal_in[gp][GlobalCoord::x];
        ny_in = intface.surface_normal_in[gp][GlobalCoord::y];

        u_in = boundary_in.q_at_gp[gp][SWE::Variables::qx] / boundary_in.aux_at_gp[gp][SWE::Auxiliaries::h];
        v_in = boundary_in.q_at_gp[gp][SWE::Variables::qy] / boundary_in.aux_at_gp[gp][SWE::Auxiliaries::h];

        uuh_in = u_in * boundary_in.q_at_gp[gp][SWE::Variables::qx];
        vvh_in = v_in * boundary_in.q_at_gp[gp][SWE::Variables::qy];
        uvh_in = u_in * boundary_in.q_at_gp[gp][SWE::Variables::qy];
        pe_in  = Global::g *
                (0.5 * std::pow(boundary_in.q_at_gp[gp][SWE::Variables::ze], 2) +
                 boundary_in.q_at_gp[gp][SWE::Variables::ze] * boundary_in.aux_at_gp[gp][SWE::Auxiliaries::bath]);

        // Fn terms
        boundary_in.F_hat_at_gp[gp][SWE::Variables::ze] =
            boundary_in.q_at_gp[gp][SWE::Variables::qx] * nx_in + boundary_in.q_at_gp[gp][SWE::Variables::qy] * ny_in;
        boundary_in.F_hat_at_gp[gp][SWE::Variables::qx] = (uuh_in + pe_in) * nx_in + uvh_in * ny_in;
        boundary_in.F_hat_at_gp[gp][SWE::Variables::qy] = uvh_in * nx_in + (vvh_in + pe_in) * ny_in;

        // dFn/dq terms
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qx) = nx_in;
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qy) = ny_in;

        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::ze) =
            (-u_in * u_in + Global::g * boundary_in.aux_at_gp[gp][SWE::Auxiliaries::h]) * nx_in - u_in * v_in * ny_in;
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qx) = 2 * u_in * nx_in + v_in * ny_in;
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qy) = u_in * ny_in;

        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::ze) =
            -u_in * v_in * nx_in + (-v_in * v_in + Global::g * boundary_in.aux_at_gp[gp][SWE::Auxiliaries::h]) * ny_in;
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qx) = v_in * nx_in;
        boundary_in.dF_hat_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qy) = u_in * nx_in + 2 * v_in * ny_in;

        /* EX State */

        nx_ex = intface.surface_normal_ex[gp_ex][GlobalCoord::x];
        ny_ex = intface.surface_normal_ex[gp_ex][GlobalCoord::y];

        u_ex = boundary_ex.q_at_gp[gp_ex][SWE::Variables::qx] / boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::h];
        v_ex = boundary_ex.q_at_gp[gp_ex][SWE::Variables::qy] / boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::h];

        uuh_ex = u_ex * boundary_ex.q_at_gp[gp_ex][SWE::Variables::qx];
        vvh_ex = v_ex * boundary_ex.q_at_gp[gp_ex][SWE::Variables::qy];
        uvh_ex = u_ex * boundary_ex.q_at_gp[gp_ex][SWE::Variables::qy];
        pe_ex  = Global::g *
                (0.5 * std::pow(boundary_ex.q_at_gp[gp_ex][SWE::Variables::ze], 2) +
                 boundary_ex.q_at_gp[gp_ex][SWE::Variables::ze] * boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::bath]);

        // Fn terms
        boundary_ex.F_hat_at_gp[gp_ex][SWE::Variables::ze] = boundary_ex.q_at_gp[gp_ex][SWE::Variables::qx] * nx_ex +
                                                             boundary_ex.q_at_gp[gp_ex][SWE::Variables::qy] * ny_ex;
        boundary_ex.F_hat_at_gp[gp_ex][SWE::Variables::qx] = (uuh_ex + pe_ex) * nx_ex + uvh_ex * ny_ex;
        boundary_ex.F_hat_at_gp[gp_ex][SWE::Variables::qy] = uvh_ex * nx_ex + (vvh_ex + pe_ex) * ny_ex;

        // dFn/dq terms
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::ze, SWE::Variables::qx) = nx_ex;
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::ze, SWE::Variables::qy) = ny_ex;

        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::qx, SWE::Variables::ze) =
            (-u_ex * u_ex + Global::g * boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::h]) * nx_ex -
            u_ex * v_ex * ny_ex;
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::qx, SWE::Variables::qx) = 2 * u_ex * nx_ex + v_ex * ny_ex;
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::qx, SWE::Variables::qy) = u_ex * ny_ex;

        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::qy, SWE::Variables::ze) =
            -u_ex * v_ex * nx_ex +
            (-v_ex * v_ex + Global::g * boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::h]) * ny_ex;
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::qy, SWE::Variables::qx) = v_ex * nx_ex;
        boundary_ex.dF_hat_dq_at_gp[gp_ex](SWE::Variables::qy, SWE::Variables::qy) = u_ex * nx_ex + 2 * v_ex * ny_ex;
    }
}
}
}

#endif
