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

    DynVector<double> q_ex(SWE::n_variables);
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        this->exchanger.GetEX(gp, q_ex);

        column(boundary.F_hat_at_gp, gp) = LLF_flux(
            Global::g, column(boundary.q_at_gp, gp), q_ex, column(boundary.aux_at_gp, gp), dbound.surface_normal[gp]);
    }

    // compute net volume flux out of IN/EX elements
    double net_volume_flux_in = 0;

    net_volume_flux_in = dbound.Integration(boundary.F_hat_at_gp)[SWE::Variables::ze];

    if (net_volume_flux_in > 0) {
        if (!wet_in) {  // water flowing from dry IN element
            // Zero flux on IN element side
            boundary.F_hat_at_gp = 0.0;

            net_volume_flux_in = 0;
        } else if (!wet_ex) {  // water flowing to dry EX element
            net_volume_flux_in = dbound.Integration(boundary.F_hat_at_gp)[SWE::Variables::ze];
        }
    } else if (net_volume_flux_in < 0) {
        if (!wet_ex) {  // water flowing from dry EX element
            // Reflective Boundary on IN element side
            BC::Land land_boundary;

            land_boundary.ComputeFlux(
                stepper, dbound.surface_normal, boundary.q_at_gp, boundary.aux_at_gp, boundary.F_hat_at_gp);

            net_volume_flux_in = 0;
        } else if (!wet_in) {  // water flowing to dry IN element
            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                this->exchanger.GetEX(gp, q_ex);

                column(boundary.F_hat_at_gp, gp) = LLF_flux(
                    0.0, column(boundary.q_at_gp, gp), q_ex, column(boundary.aux_at_gp, gp), dbound.surface_normal[gp]);
            }

            net_volume_flux_in = dbound.Integration(boundary.F_hat_at_gp)[SWE::Variables::ze];
        }
    }
}
}
}
}

#endif