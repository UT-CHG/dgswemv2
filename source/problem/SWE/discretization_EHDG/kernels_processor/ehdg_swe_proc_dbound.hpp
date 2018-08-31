#ifndef EHDG_SWE_PROC_DBOUND_HPP
#define EHDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename DistributedBoundaryType>
void Problem::global_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    boundary.q_at_gp = dbound.ComputeUgp(state.q);

    row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
        row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary state */
    auto nx = row(dbound.surface_normal, GlobalCoord::x);
    auto ny = row(dbound.surface_normal, GlobalCoord::y);

    auto u = cwise_division(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
    auto v = cwise_division(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh = cwise_multiplication(u, row(boundary.q_at_gp, SWE::Variables::qx));
    auto vvh = cwise_multiplication(v, row(boundary.q_at_gp, SWE::Variables::qy));
    auto uvh = cwise_multiplication(u, row(boundary.q_at_gp, SWE::Variables::qy));
    auto pe  = Global::g * (0.5 * cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::ze),
                                                      row(boundary.q_at_gp, SWE::Variables::ze)) +
                           cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::ze),
                                                row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary.Fn_at_gp, SWE::Variables::ze) = cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::qx), nx) +
                                                 cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::qy), ny);
    row(boundary.Fn_at_gp, SWE::Variables::qx) = cwise_multiplication(uuh + pe, nx) + cwise_multiplication(uvh, ny);
    row(boundary.Fn_at_gp, SWE::Variables::qy) = cwise_multiplication(uvh, nx) + cwise_multiplication(vvh + pe, ny);

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
    dbound.boundary_condition.exchanger.SetToSendBuffer(SWE::CommTypes::processor, message);
}

template <typename DistributedBoundaryType>
void Problem::local_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    // now compute contributions to the righthand side
    state.rhs -= dbound.IntegrationPhi(boundary.F_hat_at_gp);
}
}
}

#endif
