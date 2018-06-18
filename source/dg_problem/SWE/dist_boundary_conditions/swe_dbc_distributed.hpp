#ifndef SWE_DBC_DISTRIBUTED_HPP
#define SWE_DBC_DISTRIBUTED_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/rkdg_simulation/rkdg_stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"
#include "swe_dbc_db_data_exch.hpp"

namespace SWE {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound) {} /*nothing to initialize*/

    template <typename DistributedBoundaryType>
    void ComputeFlux(const RKDGStepper& stepper, DistributedBoundaryType& dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void Distributed::ComputeFlux(const RKDGStepper& stepper, DistributedBoundaryType& dbound) {
    bool wet_in = dbound.data.wet_dry_state.wet;
    bool wet_ex;
    this->exchanger.GetWetDryEX(wet_ex);

    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];
    auto& sp_at_gp = dbound.data.spherical_projection.sp_at_gp_boundary[dbound.bound_id];

    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        this->exchanger.GetEX(gp, ze_ex, qx_ex, qy_ex);

        LLF_flux(Global::g,
                 boundary.ze_at_gp[gp],
                 ze_ex,
                 boundary.qx_at_gp[gp],
                 qx_ex,
                 boundary.qy_at_gp[gp],
                 qy_ex,
                 boundary.bath_at_gp[gp],
                 sp_at_gp[gp],
                 dbound.surface_normal[gp],
                 boundary.ze_numerical_flux_at_gp[gp],
                 boundary.qx_numerical_flux_at_gp[gp],
                 boundary.qy_numerical_flux_at_gp[gp]);
    }

    // compute net volume flux out of IN/EX elements
    double net_volume_flux_in = 0;

    net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);

    if (net_volume_flux_in > 0) {
        if (!wet_in) {  // water flowing from dry IN element
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
            SWE::BC::Land land_boundary;

            land_boundary.ComputeFlux(stepper,
                                      dbound.surface_normal,
                                      sp_at_gp,
                                      boundary.bath_at_gp,
                                      boundary.ze_at_gp,
                                      boundary.qx_at_gp,
                                      boundary.qy_at_gp,
                                      boundary.ze_numerical_flux_at_gp,
                                      boundary.qx_numerical_flux_at_gp,
                                      boundary.qy_numerical_flux_at_gp);

            net_volume_flux_in = 0;
        } else if (!wet_in) {  // water flowing to dry IN element
            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                this->exchanger.GetEX(gp, ze_ex, qx_ex, qy_ex);

                LLF_flux(0.0,
                         boundary.ze_at_gp[gp],
                         ze_ex,
                         boundary.qx_at_gp[gp],
                         qx_ex,
                         boundary.qy_at_gp[gp],
                         qy_ex,
                         boundary.bath_at_gp[gp],
                         sp_at_gp[gp],
                         dbound.surface_normal[gp],
                         boundary.ze_numerical_flux_at_gp[gp],
                         boundary.qx_numerical_flux_at_gp[gp],
                         boundary.qy_numerical_flux_at_gp[gp]);
            }

            net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);
        }
    }
}
}
}

#endif