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
    auto pe  = Global::g * (0.5 * pow(row(boundary.q_at_gp, SWE::Variables::ze), 2.0) +
                           cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::ze),
                                                row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary.Fn_at_gp, SWE::Variables::ze) = cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::qx), nx) +
                                                 cwise_multiplication(row(boundary.q_at_gp, SWE::Variables::qy), ny);
    row(boundary.Fn_at_gp, SWE::Variables::qx) = cwise_multiplication(uuh + pe, nx) + cwise_multiplication(uvh, ny);
    row(boundary.Fn_at_gp, SWE::Variables::qy) = cwise_multiplication(uvh, nx) + cwise_multiplication(vvh + pe, ny);

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
