#ifndef EHDG_SWE_PROC_DBOUND_HPP
#define EHDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename StepperType, typename DistributedBoundaryType>
void Problem::global_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    boundary.q_at_gp = dbound.ComputeUgp(state.q);

    row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
        row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary state */
    auto nx = row(dbound.surface_normal, GlobalCoord::x);
    auto ny = row(dbound.surface_normal, GlobalCoord::y);

    auto u = vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(boundary.q_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(boundary.q_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(boundary.q_at_gp, SWE::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(boundary.q_at_gp, SWE::Variables::ze)) +
               vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary.Fn_at_gp, SWE::Variables::ze) = vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qx), nx) +
                                                 vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qy), ny);
    row(boundary.Fn_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
    row(boundary.Fn_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

    // Construct message to exterior state
    std::vector<double> message;

    message.reserve(2 * SWE::n_variables * dbound.data.get_ngp_boundary(dbound.bound_id));

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        for (uint var = 0; var < SWE::n_variables; ++var) {
            message.push_back(boundary.q_at_gp(var, gp));
        }
    }

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        for (uint var = 0; var < SWE::n_variables; ++var) {
            message.push_back(boundary.Fn_at_gp(var, gp));
        }
    }

    // Set message to send buffer
    dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::bound_state, message);
}

template <typename StepperType, typename DistributedBoundaryType>
void Problem::local_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    // now compute contributions to the righthand side
    state.rhs -= dbound.IntegrationPhi(boundary.F_hat_at_gp);
}
}
}

#endif
