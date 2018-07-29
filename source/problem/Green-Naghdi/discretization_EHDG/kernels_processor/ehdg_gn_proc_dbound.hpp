#ifndef EHDG_GN_PROC_DBOUND_HPP
#define EHDG_GN_PROC_DBOUND_HPP

namespace GN {
namespace EHDG {
template <typename DistributedBoundaryType>
void Problem::global_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    boundary.q_at_gp = dbound.ComputeUgp(state.q);

    row(boundary.aux_at_gp, GN::Auxiliaries::h) =
        row(boundary.q_at_gp, GN::Variables::ze) + row(boundary.aux_at_gp, GN::Auxiliaries::bath);

    /* Compute fluxes at boundary state */
    auto nx = row(dbound.surface_normal, GlobalCoord::x);
    auto ny = row(dbound.surface_normal, GlobalCoord::y);

    auto u = cwise_division(row(boundary.q_at_gp, GN::Variables::qx), row(boundary.aux_at_gp, GN::Auxiliaries::h));
    auto v = cwise_division(row(boundary.q_at_gp, GN::Variables::qy), row(boundary.aux_at_gp, GN::Auxiliaries::h));

    auto uuh = cwise_multiplication(u, row(boundary.q_at_gp, GN::Variables::qx));
    auto vvh = cwise_multiplication(v, row(boundary.q_at_gp, GN::Variables::qy));
    auto uvh = cwise_multiplication(u, row(boundary.q_at_gp, GN::Variables::qy));
    auto pe  = Global::g * (0.5 * cwise_multiplication(row(boundary.q_at_gp, GN::Variables::ze),
                                                      row(boundary.q_at_gp, GN::Variables::ze)) +
                           cwise_multiplication(row(boundary.q_at_gp, GN::Variables::ze),
                                                row(boundary.aux_at_gp, GN::Auxiliaries::bath)));

    row(boundary.Fn_at_gp, GN::Variables::ze) = cwise_multiplication(row(boundary.q_at_gp, GN::Variables::qx), nx) +
                                                cwise_multiplication(row(boundary.q_at_gp, GN::Variables::qy), ny);
    row(boundary.Fn_at_gp, GN::Variables::qx) = cwise_multiplication(uuh + pe, nx) + cwise_multiplication(uvh, ny);
    row(boundary.Fn_at_gp, GN::Variables::qy) = cwise_multiplication(uvh, nx) + cwise_multiplication(vvh + pe, ny);

    dbound.boundary_condition.exchanger.SetEX(boundary.q_at_gp, boundary.Fn_at_gp);
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
