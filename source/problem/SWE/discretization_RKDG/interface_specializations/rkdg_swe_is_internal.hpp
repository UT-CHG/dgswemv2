#ifndef RKDG_SWE_IS_INTERNAL_HPP
#define RKDG_SWE_IS_INTERNAL_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace ISP {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/

    template <typename InterfaceType>
    void ComputeFlux(const RKStepper& stepper, InterfaceType& intface);
};

template <typename InterfaceType>
void Internal::ComputeFlux(const RKStepper& stepper, InterfaceType& intface) {
    bool wet_in = intface.data_in.wet_dry_state.wet;
    bool wet_ex = intface.data_ex.wet_dry_state.wet;

    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    // assemble numerical fluxes
    uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
    uint gp_ex;
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        gp_ex = ngp - gp - 1;

        column(boundary_in.F_hat_at_gp, gp) = LLF_flux(Global::g,
                                                       column(boundary_in.q_at_gp, gp),
                                                       column(boundary_ex.q_at_gp, gp_ex),
                                                       column(boundary_in.aux_at_gp, gp),
                                                       column(intface.surface_normal_in, gp));

        column(boundary_ex.F_hat_at_gp, gp_ex) = -column(boundary_in.F_hat_at_gp, gp);
    }

    // compute net volume flux out of IN/EX elements
    double net_volume_flux_in = 0;
    // double net_volume_flux_ex = 0;

    net_volume_flux_in = intface.IntegrationIN(boundary_in.F_hat_at_gp)[SWE::Variables::ze];
    // net_volume_flux_ex = -net_volume_flux_in;

    if (net_volume_flux_in > 0) {
        if (!wet_in) {  // water flowing from dry IN element
            // Zero flux on IN element side
            set_constant(boundary_in.F_hat_at_gp, 0.0);

            // Reflective Boundary on EX element side
            BC::Land land_boundary;
            land_boundary.Initialize(intface.data_ex.get_ngp_boundary(intface.bound_id_ex));

            land_boundary.ComputeFlux(stepper,
                                      intface.surface_normal_ex,
                                      boundary_ex.q_at_gp,
                                      boundary_ex.aux_at_gp,
                                      boundary_ex.F_hat_at_gp);

            net_volume_flux_in = 0;
            // net_volume_flux_ex = 0;
        } else if (!wet_ex) {  // water flowing to dry EX element
            for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                gp_ex = ngp - gp - 1;

                column(boundary_ex.F_hat_at_gp, gp_ex) = LLF_flux(0.0,
                                                                  column(boundary_ex.q_at_gp, gp_ex),
                                                                  column(boundary_in.q_at_gp, gp),
                                                                  column(boundary_ex.aux_at_gp, gp_ex),
                                                                  column(intface.surface_normal_ex, gp_ex));
            }

            net_volume_flux_in = intface.IntegrationIN(boundary_in.F_hat_at_gp)[SWE::Variables::ze];
            // net_volume_flux_ex = intface.IntegrationEX(boundary_ex.ze_numerical_flux_at_gp);
        }
    } else if (net_volume_flux_in < 0) {
        if (!wet_ex) {  // water flowing from dry EX element
            // Zero flux on EX element side
            set_constant(boundary_ex.F_hat_at_gp, 0.0);

            // Reflective Boundary on IN element side
            BC::Land land_boundary;
            land_boundary.Initialize(intface.data_in.get_ngp_boundary(intface.bound_id_in));

            land_boundary.ComputeFlux(stepper,
                                      intface.surface_normal_in,
                                      boundary_in.q_at_gp,
                                      boundary_in.aux_at_gp,
                                      boundary_in.F_hat_at_gp);

            net_volume_flux_in = 0;
            // net_volume_flux_ex = 0;
        } else if (!wet_in) {  // water flowing to dry IN element
            for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                gp_ex = ngp - gp - 1;

                column(boundary_in.F_hat_at_gp, gp) = LLF_flux(0.0,
                                                               column(boundary_in.q_at_gp, gp),
                                                               column(boundary_ex.q_at_gp, gp_ex),
                                                               column(boundary_in.aux_at_gp, gp),
                                                               column(intface.surface_normal_in, gp));
            }

            net_volume_flux_in = intface.IntegrationIN(boundary_in.F_hat_at_gp)[SWE::Variables::ze];
            // net_volume_flux_ex = intface.IntegrationEX(boundary_ex.ze_numerical_flux_at_gp);
        }
    }
}
}
}
}

#endif