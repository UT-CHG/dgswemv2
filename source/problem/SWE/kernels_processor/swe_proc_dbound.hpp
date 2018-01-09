#ifndef SWE_PROC_DBOUND_HPP
#define SWE_PROC_DBOUND_HPP

#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
template <typename DistributedBoundaryType>
void Problem::distributed_boundary_send_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.get_stage();

    auto& state = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.ze, boundary.ze_at_gp);
    dbound.ComputeUgp(state.qx, boundary.qx_at_gp);
    dbound.ComputeUgp(state.qy, boundary.qy_at_gp);

    dbound.boundary_condition.SetEX(boundary.ze_at_gp, boundary.qx_at_gp, boundary.qy_at_gp);

    dbound.boundary_condition.SetWetDryEX(dbound.data.wet_dry_state.wet);
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    auto& wd_state_in = dbound.data.wet_dry_state;
    
    bool wet_ex;
    dbound.boundary_condition.GetWetDryEX(wet_ex);

    if (wd_state_in.wet || wet_ex) {
        const uint stage = stepper.get_stage();
        const double dt = stepper.get_dt();

        auto& state = dbound.data.state[stage];
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        double ze_ex, qx_ex, qy_ex;
        for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
            dbound.boundary_condition.GetEX(stepper,
                                            gp,
                                            dbound.surface_normal,
                                            boundary.ze_at_gp,
                                            boundary.qx_at_gp,
                                            boundary.qy_at_gp,
                                            ze_ex,
                                            qx_ex,
                                            qy_ex);

            LLF_flux(boundary.ze_at_gp[gp],
                     ze_ex,
                     boundary.qx_at_gp[gp],
                     qx_ex,
                     boundary.qy_at_gp[gp],
                     qy_ex,
                     boundary.bath_at_gp[gp],
                     dbound.surface_normal[gp],
                     boundary.ze_numerical_flux_at_gp[gp],
                     boundary.qx_numerical_flux_at_gp[gp],
                     boundary.qy_numerical_flux_at_gp[gp]);
        }

        // compute net volume flux out of IN/EX elements
        double net_volume_flux_in = 0;

        net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);

        if (net_volume_flux_in > 0) {
            if (!wd_state_in.wet) {  // water flowing from dry IN element
                // Zero flux on IN element side
                std::fill(boundary.ze_numerical_flux_at_gp.begin(), boundary.ze_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary.qx_numerical_flux_at_gp.begin(), boundary.qx_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary.qy_numerical_flux_at_gp.begin(), boundary.qy_numerical_flux_at_gp.end(), 0.0);

                net_volume_flux_in = 0;
            } else if (!wet_ex) {  // water flowing to dry EX element
                net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);
            }
        } else if (net_volume_flux_in < 0) {
            if (!wet_ex) {  // water flowing from dry EX element
                // Reflective Boundary on IN element side
                SWE::Land land_boundary;

                for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                    land_boundary.GetEX(stepper,
                                        gp,
                                        dbound.surface_normal,
                                        boundary.ze_at_gp,
                                        boundary.qx_at_gp,
                                        boundary.qy_at_gp,
                                        ze_ex,
                                        qx_ex,
                                        qy_ex);

                    LLF_flux(boundary.ze_at_gp[gp],
                             ze_ex,
                             boundary.qx_at_gp[gp],
                             qx_ex,
                             boundary.qy_at_gp[gp],
                             qy_ex,
                             boundary.bath_at_gp[gp],
                             dbound.surface_normal[gp],
                             boundary.ze_numerical_flux_at_gp[gp],
                             boundary.qx_numerical_flux_at_gp[gp],
                             boundary.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = 0;
            } else if (!wd_state_in.wet) {  // water flowing to dry IN element
                for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                    dbound.boundary_condition.GetEX(stepper,
                                            gp,
                                            dbound.surface_normal,
                                            boundary.ze_at_gp,
                                            boundary.qx_at_gp,
                                            boundary.qy_at_gp,
                                            ze_ex,
                                            qx_ex,
                                            qy_ex);

                    LLF_flux_zero_g(boundary.ze_at_gp[gp],
                                    ze_ex,
                                    boundary.qx_at_gp[gp],
                                    qx_ex,
                                    boundary.qy_at_gp[gp],
                                    qy_ex,
                                    boundary.bath_at_gp[gp],
                                    dbound.surface_normal[gp],
                                    boundary.ze_numerical_flux_at_gp[gp],
                                    boundary.qx_numerical_flux_at_gp[gp],
                                    boundary.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);
            }
        }

        // compute water level reduction
        wd_state_in.water_volume -= net_volume_flux_in * dt;

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < dbound.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] -= dbound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
            state.rhs_qx[dof] -= dbound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
            state.rhs_qy[dof] -= dbound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
        }
    }
}
}

#endif
