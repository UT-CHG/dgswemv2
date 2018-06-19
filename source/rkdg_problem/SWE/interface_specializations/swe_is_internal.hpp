#ifndef SWE_IS_INTERNAL_HPP
#define SWE_IS_INTERNAL_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/rkdg_simulation/rk_stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace IS {
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

    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& sp_at_gp_in = intface.data_in.spherical_projection.sp_at_gp_boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
    auto& sp_at_gp_ex = intface.data_ex.spherical_projection.sp_at_gp_boundary[intface.bound_id_ex];

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
        if (!wet_in) {  // water flowing from dry IN element
            // Zero flux on IN element side
            std::fill(boundary_in.ze_numerical_flux_at_gp.begin(), boundary_in.ze_numerical_flux_at_gp.end(), 0.0);
            std::fill(boundary_in.qx_numerical_flux_at_gp.begin(), boundary_in.qx_numerical_flux_at_gp.end(), 0.0);
            std::fill(boundary_in.qy_numerical_flux_at_gp.begin(), boundary_in.qy_numerical_flux_at_gp.end(), 0.0);

            // Reflective Boundary on EX element side
            SWE::BC::Land land_boundary;

            land_boundary.ComputeFlux(stepper,
                                      intface.surface_normal_ex,
                                      sp_at_gp_ex,
                                      boundary_ex.bath_at_gp,
                                      boundary_ex.ze_at_gp,
                                      boundary_ex.qx_at_gp,
                                      boundary_ex.qy_at_gp,
                                      boundary_ex.ze_numerical_flux_at_gp,
                                      boundary_ex.qx_numerical_flux_at_gp,
                                      boundary_ex.qy_numerical_flux_at_gp);

            net_volume_flux_in = 0;
            net_volume_flux_ex = 0;
        } else if (!wet_ex) {  // water flowing to dry EX element
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
        if (!wet_ex) {  // water flowing from dry EX element
            // Zero flux on EX element side
            std::fill(boundary_ex.ze_numerical_flux_at_gp.begin(), boundary_ex.ze_numerical_flux_at_gp.end(), 0.0);
            std::fill(boundary_ex.qx_numerical_flux_at_gp.begin(), boundary_ex.qx_numerical_flux_at_gp.end(), 0.0);
            std::fill(boundary_ex.qy_numerical_flux_at_gp.begin(), boundary_ex.qy_numerical_flux_at_gp.end(), 0.0);

            // Reflective Boundary on IN element side
            SWE::BC::Land land_boundary;

            land_boundary.ComputeFlux(stepper,
                                      intface.surface_normal_in,
                                      sp_at_gp_in,
                                      boundary_in.bath_at_gp,
                                      boundary_in.ze_at_gp,
                                      boundary_in.qx_at_gp,
                                      boundary_in.qy_at_gp,
                                      boundary_in.ze_numerical_flux_at_gp,
                                      boundary_in.qx_numerical_flux_at_gp,
                                      boundary_in.qy_numerical_flux_at_gp);

            net_volume_flux_in = 0;
            net_volume_flux_ex = 0;
        } else if (!wet_in) {  // water flowing to dry IN element
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
}
}
}

#endif