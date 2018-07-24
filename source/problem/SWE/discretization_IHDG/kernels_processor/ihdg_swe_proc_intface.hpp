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
    auto pe_in  = Global::g * (0.5 * cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::ze),
                                                         row(boundary_in.q_at_gp, SWE::Variables::ze)) +
                              cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::ze),
                                                   row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath)));

    // Fn terms
    row(boundary_in.F_hat_at_gp, SWE::Variables::ze) =
        cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::qx), nx_in) +
        cwise_multiplication(row(boundary_in.q_at_gp, SWE::Variables::qy), ny_in);
    row(boundary_in.F_hat_at_gp, SWE::Variables::qx) =
        cwise_multiplication(uuh_in + pe_in, nx_in) + cwise_multiplication(uvh_in, ny_in);
    row(boundary_in.F_hat_at_gp, SWE::Variables::qy) =
        cwise_multiplication(uvh_in, nx_in) + cwise_multiplication(vvh_in + pe_in, ny_in);

    // dFn/dq terms
    set_constant(row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::ze_ze), 0.0);
    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::ze_qx) = nx_in;
    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::ze_qy) = ny_in;

    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::qx_ze) =
        cwise_multiplication(
            -cwise_multiplication(u_in, u_in) + Global::g * row(boundary_in.aux_at_gp, SWE::Auxiliaries::h), nx_in) -
        cwise_multiplication(cwise_multiplication(u_in, v_in), ny_in);
    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::qx_qx) =
        2.0 * cwise_multiplication(u_in, nx_in) + cwise_multiplication(v_in, ny_in);
    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::qx_qy) = cwise_multiplication(u_in, ny_in);

    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::qy_ze) =
        -cwise_multiplication(cwise_multiplication(u_in, v_in), nx_in) +
        cwise_multiplication(
            -cwise_multiplication(v_in, v_in) + Global::g * row(boundary_in.aux_at_gp, SWE::Auxiliaries::h), ny_in);
    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::qy_qx) = cwise_multiplication(v_in, nx_in);
    row(boundary_in.dF_hat_dq_at_gp, JacobianVariables::qy_qy) =
        cwise_multiplication(u_in, nx_in) + 2.0 * cwise_multiplication(v_in, ny_in);

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
    auto pe_ex  = Global::g * (0.5 * cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::ze),
                                                         row(boundary_ex.q_at_gp, SWE::Variables::ze)) +
                              cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::ze),
                                                   row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath)));

    // Fn terms
    row(boundary_ex.F_hat_at_gp, SWE::Variables::ze) =
        cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::qx), nx_ex) +
        cwise_multiplication(row(boundary_ex.q_at_gp, SWE::Variables::qy), ny_ex);
    row(boundary_ex.F_hat_at_gp, SWE::Variables::qx) =
        cwise_multiplication(uuh_ex + pe_ex, nx_ex) + cwise_multiplication(uvh_ex, ny_ex);
    row(boundary_ex.F_hat_at_gp, SWE::Variables::qy) =
        cwise_multiplication(uvh_ex, nx_ex) + cwise_multiplication(vvh_ex + pe_ex, ny_ex);

    // dFn/dq terms
    set_constant(row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::ze_ze), 0.0);
    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::ze_qx) = nx_ex;
    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::ze_qy) = ny_ex;

    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::qx_ze) =
        cwise_multiplication(
            -cwise_multiplication(u_ex, u_ex) + Global::g * row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h), nx_ex) -
        cwise_multiplication(cwise_multiplication(u_ex, v_ex), ny_ex);
    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::qx_qx) =
        2.0 * cwise_multiplication(u_ex, nx_ex) + cwise_multiplication(v_ex, ny_ex);
    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::qx_qy) = cwise_multiplication(u_ex, ny_ex);

    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::qy_ze) =
        -cwise_multiplication(cwise_multiplication(u_ex, v_ex), nx_ex) +
        cwise_multiplication(
            -cwise_multiplication(v_ex, v_ex) + Global::g * row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h), ny_ex);
    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::qy_qx) = cwise_multiplication(v_ex, nx_ex);
    row(boundary_ex.dF_hat_dq_at_gp, JacobianVariables::qy_qy) =
        cwise_multiplication(u_ex, nx_ex) + 2.0 * cwise_multiplication(v_ex, ny_ex);
}
}
}

#endif
