#ifndef RKDG_SWE_DBC_DISTRIBUTED_HPP
#define RKDG_SWE_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"
#include "rkdg_swe_dbc_db_data_exch.hpp"

namespace SWE {
namespace RKDG {
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
    void ComputeFlux(const RKStepper& stepper, DistributedBoundaryType& dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void Distributed::ComputeFlux(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    bool wet_in = dbound.data.wet_dry_state.wet;
    bool wet_ex;
    this->exchanger.GetWetDryEX(wet_ex);

    auto& boundary = dbound.data.boundary[dbound.bound_id];
    auto& sp_at_gp = dbound.data.spherical_projection.sp_at_gp_boundary[dbound.bound_id];

    Vector<double, SWE::n_variables> q_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        this->exchanger.GetEX(gp, q_ex);

        LLF_flux(Global::g,
                 boundary.q_at_gp[gp],
                 q_ex,
                 boundary.bath_at_gp[gp],
                 sp_at_gp[gp],
                 dbound.surface_normal[gp],
                 boundary.F_hat_at_gp[gp]);
    }

    // compute net volume flux out of IN/EX elements
    double net_volume_flux_in = 0;

    net_volume_flux_in = dbound.Integration(boundary.F_hat_at_gp)[SWE::Variables::ze];

    if (net_volume_flux_in > 0) {
        if (!wet_in) {  // water flowing from dry IN element
            // Zero flux on IN element side
            std::fill(boundary.F_hat_at_gp.begin(), boundary.F_hat_at_gp.end(), 0.0);

            net_volume_flux_in = 0;
        } else if (!wet_ex) {  // water flowing to dry EX element
            net_volume_flux_in = dbound.Integration(boundary.F_hat_at_gp)[SWE::Variables::ze];
        }
    } else if (net_volume_flux_in < 0) {
        if (!wet_ex) {  // water flowing from dry EX element
            // Reflective Boundary on IN element side
            BC::Land land_boundary;

            land_boundary.ComputeFlux(
                stepper, dbound.surface_normal, sp_at_gp, boundary.bath_at_gp, boundary.q_at_gp, boundary.F_hat_at_gp);

            net_volume_flux_in = 0;
        } else if (!wet_in) {  // water flowing to dry IN element
            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                this->exchanger.GetEX(gp, q_ex);

                LLF_flux(0.0,
                         boundary.q_at_gp[gp],
                         q_ex,
                         boundary.bath_at_gp[gp],
                         sp_at_gp[gp],
                         dbound.surface_normal[gp],
                         boundary.F_hat_at_gp[gp]);
            }

            net_volume_flux_in = dbound.Integration(boundary.F_hat_at_gp)[SWE::Variables::ze];
        }
    }
}
}
}
}

#endif