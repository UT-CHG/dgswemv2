#ifndef EHDG_SWE_PROC_INTFACE_HPP
#define EHDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace EHDG {
template <typename InterfaceType>
void Problem::global_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
    boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);

    uint gp_ex;
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        gp_ex = intface.data_in.get_ngp_boundary(intface.bound_id_in) - gp - 1;

        boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp) =
            boundary_in.q_at_gp(SWE::Variables::ze, gp) + boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp);
        boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex) =
            boundary_ex.q_at_gp(SWE::Variables::ze, gp_ex) + boundary_ex.aux_at_gp(SWE::Auxiliaries::bath, gp_ex);
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

        nx_in = intface.surface_normal_in(GlobalCoord::x, gp);
        ny_in = intface.surface_normal_in(GlobalCoord::y, gp);

        u_in = boundary_in.q_at_gp(SWE::Variables::qx, gp) / boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp);
        v_in = boundary_in.q_at_gp(SWE::Variables::qy, gp) / boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp);

        uuh_in = u_in * boundary_in.q_at_gp(SWE::Variables::qx, gp);
        vvh_in = v_in * boundary_in.q_at_gp(SWE::Variables::qy, gp);
        uvh_in = u_in * boundary_in.q_at_gp(SWE::Variables::qy, gp);
        pe_in  = Global::g *
                (0.5 * std::pow(boundary_in.q_at_gp(SWE::Variables::ze, gp), 2) +
                 boundary_in.q_at_gp(SWE::Variables::ze, gp) * boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp));

        boundary_in.Fn_at_gp(SWE::Variables::ze, gp) =
            boundary_in.q_at_gp(SWE::Variables::qx, gp) * nx_in + boundary_in.q_at_gp(SWE::Variables::qy, gp) * ny_in;
        boundary_in.Fn_at_gp(SWE::Variables::qx, gp) = (uuh_in + pe_in) * nx_in + uvh_in * ny_in;
        boundary_in.Fn_at_gp(SWE::Variables::qy, gp) = uvh_in * nx_in + (vvh_in + pe_in) * ny_in;

        /* EX State */

        nx_ex = intface.surface_normal_ex(GlobalCoord::x, gp_ex);
        ny_ex = intface.surface_normal_ex(GlobalCoord::y, gp_ex);

        u_ex = boundary_ex.q_at_gp(SWE::Variables::qx, gp_ex) / boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex);
        v_ex = boundary_ex.q_at_gp(SWE::Variables::qy, gp_ex) / boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex);

        uuh_ex = u_ex * boundary_ex.q_at_gp(SWE::Variables::qx, gp_ex);
        vvh_ex = v_ex * boundary_ex.q_at_gp(SWE::Variables::qy, gp_ex);
        uvh_ex = u_ex * boundary_ex.q_at_gp(SWE::Variables::qy, gp_ex);
        pe_ex  = Global::g *
                (0.5 * std::pow(boundary_ex.q_at_gp(SWE::Variables::ze, gp_ex), 2) +
                 boundary_ex.q_at_gp(SWE::Variables::ze, gp_ex) * boundary_ex.aux_at_gp(SWE::Auxiliaries::bath, gp_ex));

        boundary_ex.Fn_at_gp(SWE::Variables::ze, gp_ex) = boundary_ex.q_at_gp(SWE::Variables::qx, gp_ex) * nx_ex +
                                                          boundary_ex.q_at_gp(SWE::Variables::qy, gp_ex) * ny_ex;
        boundary_ex.Fn_at_gp(SWE::Variables::qx, gp_ex) = (uuh_ex + pe_ex) * nx_ex + uvh_ex * ny_ex;
        boundary_ex.Fn_at_gp(SWE::Variables::qy, gp_ex) = uvh_ex * nx_ex + (vvh_ex + pe_ex) * ny_ex;
    }
}

template <typename InterfaceType>
void Problem::local_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
        column(state_in.rhs, dof) -= intface.IntegrationPhiIN(dof, boundary_in.F_hat_at_gp);
    }

    for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
        column(state_ex.rhs, dof) -= intface.IntegrationPhiEX(dof, boundary_ex.F_hat_at_gp);
    }
}
}
}

#endif
