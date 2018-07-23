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

    row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
        row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);

    row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
        row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary states */

    /* IN State */
    auto nx_in = row(intface.surface_normal_in, GlobalCoord::x);
    auto ny_in = row(intface.surface_normal_in, GlobalCoord::y);

    auto u_in =
        cwise_division(row(boundary_in.q_at_gp, SWE::Variables::qx), row(boundary_in.aux_at_gp, SWE::Auxiliaries::h));
    auto v_in =
        cwise_division(row(boundary_in.q_at_gp, SWE::Variables::qy), row(boundary_in.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh_in = cwise_multiplication(u_in, row(boundary_in.q_at_gp, SWE::Variables::qx));
    auto vvh_in = cwise_multiplication(v_in, row(boundary_in.q_at_gp, SWE::Variables::qy));
    auto uvh_in = cwise_multiplication(u_in, row(boundary_in.q_at_gp, SWE::Variables::qy));
    auto pe_in  = Global::g * (0.5 * pow(row(boundary_in.q_at_gp, SWE::Variables::ze), 2.0) +
                              cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::ze),
                                                   row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary_in.Fn_at_gp, SWE::Variables::ze) =
        cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::qx), nx_in) +
        cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::qy), ny_in);
    row(boundary_in.Fn_at_gp, SWE::Variables::qx) =
        cwise_multiplication(uuh_in + pe_in, nx_in) + cwise_multiplication(uvh_in, ny_in);
    row(boundary_in.Fn_at_gp, SWE::Variables::qy) =
        cwise_multiplication(uvh_in, nx_in) + cwise_multiplication(vvh_in + pe_in, ny_in);

    /* EX State */
    auto nx_ex = row(intface.surface_normal_ex, GlobalCoord::x);
    auto ny_ex = row(intface.surface_normal_ex, GlobalCoord::y);

    auto u_ex =
        cwise_division(row(boundary_ex.q_at_gp, SWE::Variables::qx), row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h));
    auto v_ex =
        cwise_division(row(boundary_ex.q_at_gp, SWE::Variables::qy), row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh_ex = cwise_multiplication(u_ex, row(boundary_ex.q_at_gp, SWE::Variables::qx));
    auto vvh_ex = cwise_multiplication(v_ex, row(boundary_ex.q_at_gp, SWE::Variables::qy));
    auto uvh_ex = cwise_multiplication(u_ex, row(boundary_ex.q_at_gp, SWE::Variables::qy));
    auto pe_ex  = Global::g * (0.5 * pow(row(boundary_ex.q_at_gp, SWE::Variables::ze), 2.0) +
                              cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::ze),
                                                   row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary_ex.Fn_at_gp, SWE::Variables::ze) =
        cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::qx), nx_ex) +
        cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::qy), ny_ex);
    row(boundary_ex.Fn_at_gp, SWE::Variables::qx) =
        cwise_multiplication(uuh_ex + pe_ex, nx_ex) + cwise_multiplication(uvh_ex, ny_ex);
    row(boundary_ex.Fn_at_gp, SWE::Variables::qy) =
        cwise_multiplication(uvh_ex, nx_ex) + cwise_multiplication(vvh_ex + pe_ex, ny_ex);
}

template <typename InterfaceType>
void Problem::local_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    // now compute contributions to the righthand side
    state_in.rhs -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp);

    state_ex.rhs -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp);
}
}
}

#endif
