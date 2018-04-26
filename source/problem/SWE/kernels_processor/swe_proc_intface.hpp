#ifndef SWE_PROC_INTFACE_HPP
#define SWE_PROC_INTFACE_HPP

#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
template <typename InterfaceType>
void Problem::interface_kernel(const Stepper& stepper, InterfaceType& intface) {
    auto& wd_state_in = intface.data_in.wet_dry_state;
    auto& wd_state_ex = intface.data_ex.wet_dry_state;

    if (wd_state_in.wet || wd_state_ex.wet) {
        const uint   stage = stepper.GetStage();
        const double dt    = stepper.GetDT();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& sp_at_gp_in = intface.data_in.spherical_projection.sp_at_gp_boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        auto& sp_at_gp_ex = intface.data_ex.spherical_projection.sp_at_gp_boundary[intface.bound_id_ex];

        intface.ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
        intface.ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
        intface.ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

        intface.ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
        intface.ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
        intface.ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);

        // assemble numerical fluxes
        uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex = 0;
        for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
            gp_ex = ngp - gp - 1;

            LLF_flux(Global::g,
                     boundary_in.ze_at_gp[gp],
                     boundary_ex.ze_at_gp[gp_ex],
                     boundary_in.qx_at_gp[gp],
                     boundary_ex.qx_at_gp[gp_ex],
                     boundary_in.qy_at_gp[gp],
                     boundary_ex.qy_at_gp[gp_ex],
                     boundary_in.bath_at_gp[gp],
                     sp_at_gp_in[gp],
                     intface.surface_normal_in[gp],
                     boundary_in.ze_numerical_flux_at_gp[gp],
                     boundary_in.qx_numerical_flux_at_gp[gp],
                     boundary_in.qy_numerical_flux_at_gp[gp]);

            boundary_ex.ze_numerical_flux_at_gp[gp_ex] = -boundary_in.ze_numerical_flux_at_gp[gp];
            boundary_ex.qx_numerical_flux_at_gp[gp_ex] = -boundary_in.qx_numerical_flux_at_gp[gp];
            boundary_ex.qy_numerical_flux_at_gp[gp_ex] = -boundary_in.qy_numerical_flux_at_gp[gp];
        }

        // compute net volume flux out of IN/EX elements
        double net_volume_flux_in = 0;
        double net_volume_flux_ex = 0;

        net_volume_flux_in = intface.IntegrationIN(boundary_in.ze_numerical_flux_at_gp);
        net_volume_flux_ex = -net_volume_flux_in;

        if (net_volume_flux_in > 0) {
            if (!wd_state_in.wet) {  // water flowing from dry IN element
                // Zero flux on IN element side
                std::fill(boundary_in.ze_numerical_flux_at_gp.begin(), boundary_in.ze_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_in.qx_numerical_flux_at_gp.begin(), boundary_in.qx_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_in.qy_numerical_flux_at_gp.begin(), boundary_in.qy_numerical_flux_at_gp.end(), 0.0);

                // Reflective Boundary on EX element side
                SWE::Land land_boundary;
                double    ze_in, qx_in, qy_in;

                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    land_boundary.GetEX(stepper,
                                        gp,
                                        intface.surface_normal_ex,
                                        boundary_ex.ze_at_gp,
                                        boundary_ex.qx_at_gp,
                                        boundary_ex.qy_at_gp,
                                        ze_in,
                                        qx_in,
                                        qy_in);

                    LLF_flux(Global::g,
                             boundary_ex.ze_at_gp[gp],
                             ze_in,
                             boundary_ex.qx_at_gp[gp],
                             qx_in,
                             boundary_ex.qy_at_gp[gp],
                             qy_in,
                             boundary_ex.bath_at_gp[gp],
                             sp_at_gp_ex[gp],
                             intface.surface_normal_ex[gp],
                             boundary_ex.ze_numerical_flux_at_gp[gp],
                             boundary_ex.qx_numerical_flux_at_gp[gp],
                             boundary_ex.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = 0;
                net_volume_flux_ex = 0;
            } else if (!wd_state_ex.wet) {  // water flowing to dry EX element
                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    gp_ex = ngp - gp - 1;

                    LLF_flux(0.0,
                             boundary_ex.ze_at_gp[gp_ex],
                             boundary_in.ze_at_gp[gp],
                             boundary_ex.qx_at_gp[gp_ex],
                             boundary_in.qx_at_gp[gp],
                             boundary_ex.qy_at_gp[gp_ex],
                             boundary_in.qy_at_gp[gp],
                             boundary_ex.bath_at_gp[gp_ex],
                             sp_at_gp_ex[gp_ex],
                             intface.surface_normal_ex[gp_ex],
                             boundary_ex.ze_numerical_flux_at_gp[gp_ex],
                             boundary_ex.qx_numerical_flux_at_gp[gp_ex],
                             boundary_ex.qy_numerical_flux_at_gp[gp_ex]);
                }

                net_volume_flux_in = intface.IntegrationIN(boundary_in.ze_numerical_flux_at_gp);
                net_volume_flux_ex = intface.IntegrationEX(boundary_ex.ze_numerical_flux_at_gp);
            }
        } else if (net_volume_flux_in < 0) {
            if (!wd_state_ex.wet) {  // water flowing from dry EX element
                // Zero flux on EX element side
                std::fill(boundary_ex.ze_numerical_flux_at_gp.begin(), boundary_ex.ze_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_ex.qx_numerical_flux_at_gp.begin(), boundary_ex.qx_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_ex.qy_numerical_flux_at_gp.begin(), boundary_ex.qy_numerical_flux_at_gp.end(), 0.0);

                // Reflective Boundary on IN element side
                SWE::Land land_boundary;
                double    ze_ex, qx_ex, qy_ex;

                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    land_boundary.GetEX(stepper,
                                        gp,
                                        intface.surface_normal_in,
                                        boundary_in.ze_at_gp,
                                        boundary_in.qx_at_gp,
                                        boundary_in.qy_at_gp,
                                        ze_ex,
                                        qx_ex,
                                        qy_ex);

                    LLF_flux(Global::g,
                             boundary_in.ze_at_gp[gp],
                             ze_ex,
                             boundary_in.qx_at_gp[gp],
                             qx_ex,
                             boundary_in.qy_at_gp[gp],
                             qy_ex,
                             boundary_in.bath_at_gp[gp],
                             sp_at_gp_in[gp],
                             intface.surface_normal_in[gp],
                             boundary_in.ze_numerical_flux_at_gp[gp],
                             boundary_in.qx_numerical_flux_at_gp[gp],
                             boundary_in.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = 0;
                net_volume_flux_ex = 0;
            } else if (!wd_state_in.wet) {  // water flowing to dry IN element
                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    gp_ex = ngp - gp - 1;

                    LLF_flux(0.0,
                             boundary_in.ze_at_gp[gp],
                             boundary_ex.ze_at_gp[gp_ex],
                             boundary_in.qx_at_gp[gp],
                             boundary_ex.qx_at_gp[gp_ex],
                             boundary_in.qy_at_gp[gp],
                             boundary_ex.qy_at_gp[gp_ex],
                             boundary_in.bath_at_gp[gp],
                             sp_at_gp_in[gp],
                             intface.surface_normal_in[gp],
                             boundary_in.ze_numerical_flux_at_gp[gp],
                             boundary_in.qx_numerical_flux_at_gp[gp],
                             boundary_in.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = intface.IntegrationIN(boundary_in.ze_numerical_flux_at_gp);
                net_volume_flux_ex = intface.IntegrationEX(boundary_ex.ze_numerical_flux_at_gp);
            }
        }

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
            state_in.rhs_ze[dof] -= intface.IntegrationPhiIN(dof, boundary_in.ze_numerical_flux_at_gp);
            state_in.rhs_qx[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qx_numerical_flux_at_gp);
            state_in.rhs_qy[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qy_numerical_flux_at_gp);
        }

        for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
            state_ex.rhs_ze[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.ze_numerical_flux_at_gp);
            state_ex.rhs_qx[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qx_numerical_flux_at_gp);
            state_ex.rhs_qy[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qy_numerical_flux_at_gp);
        }
    }
}
}

#endif
