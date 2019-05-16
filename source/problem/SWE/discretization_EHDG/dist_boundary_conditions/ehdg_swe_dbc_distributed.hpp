#ifndef EHDG_SWE_DBC_DISTRIBUTED_HPP
#define EHDG_SWE_DBC_DISTRIBUTED_HPP

#include "communication/db_data_exchanger.hpp"

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

  private:
    HybMatrix<double, SWE::n_variables> q_ex;

    BC::Land land_boundary;

  public:
    Distributed(DBDataExchanger  exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(DBDataExchanger  exchanger) : exchanger(std::move(exchanger)) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {
    uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

    this->q_ex.resize(SWE::n_variables, ngp);
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    /* Get message from ex */

    std::vector<double> message;

    uint ngp = edge_dbound.edge_data.get_ngp();

    message.resize(1 + SWE::n_variables * ngp);

    edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);

    bool wet_ex = (bool)message[0];

    uint gp_ex;
    for (uint gp = 0; gp < ngp; ++gp) {
        gp_ex = ngp - gp - 1;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            this->q_ex(var, gp_ex) = message[1 + SWE::n_variables * gp + var];
        }
    }

    bool wet_in = edge_dbound.boundary.data.wet_dry_state.wet;

    if (wet_in || wet_ex) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        // Our definition of numerical flux implies q_hat = 0.5 * (q_in + q_ex)
        edge_internal.q_hat_at_gp = (boundary.q_at_gp + this->q_ex) / 2.0;

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        /* Compute trace flux */

        SWE::get_Fn(edge_internal.q_hat_at_gp,
                    edge_internal.aux_hat_at_gp,
                    edge_dbound.boundary.surface_normal,
                    boundary.F_hat_at_gp);

        /* Add stabilization parameter terms */

        SWE::get_tau_LF(edge_internal.q_hat_at_gp,
                        edge_internal.aux_hat_at_gp,
                        edge_dbound.boundary.surface_normal,
                        edge_internal.tau);

        for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
            column(boundary.F_hat_at_gp, gp) +=
                edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
        }

        // correct numerical fluxes for wetting/drying
        for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
            if (boundary.F_hat_at_gp(Variables::ze, gp) > 1e-12) {
                if (!wet_in) {  // water flowing from dry IN element
                    // Zero flux on IN element side
                    set_constant(column(boundary.F_hat_at_gp, gp), 0.0);
                }
            } else if (boundary.F_hat_at_gp(Variables::ze, gp) < -1e-12) {
                if (!wet_ex) {  // water flowing from dry EX element
                    // Reflective Boundary on IN element side

                    this->land_boundary.ComputeFlux(column(edge_dbound.boundary.surface_normal, gp),
                                                    column(boundary.q_at_gp, gp),
                                                    column(edge_internal.aux_hat_at_gp, gp),
                                                    edge_internal.tau[gp],
                                                    column(edge_internal.q_hat_at_gp, gp),
                                                    column(boundary.F_hat_at_gp, gp));

                } else if (!wet_in) {  // water flowing to dry IN element
                    double temp_flux_ze = boundary.F_hat_at_gp(Variables::ze, gp);

                    SWE::get_Fn(0.0,
                                column(edge_internal.q_hat_at_gp, gp),
                                column(edge_internal.aux_hat_at_gp, gp),
                                column(edge_dbound.boundary.surface_normal, gp),
                                column(boundary.F_hat_at_gp, gp));

                    SWE::get_tau_LF(0.0,
                                    column(edge_internal.q_hat_at_gp, gp),
                                    column(edge_internal.aux_hat_at_gp, gp),
                                    column(edge_dbound.boundary.surface_normal, gp),
                                    edge_internal.tau[gp]);

                    column(boundary.F_hat_at_gp, gp) +=
                        edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));

                    // Only remove gravity contributions for the momentum fluxes
                    boundary.F_hat_at_gp(Variables::ze, gp) = temp_flux_ze;
                }
            }

            assert(!std::isnan(boundary.F_hat_at_gp(Variables::ze, gp)));
        }
    }
}
}
}
}

#endif