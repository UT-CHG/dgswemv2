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
    HybMatrix<double, SWE::n_auxiliaries> aux_ex;

    BC::Land land_boundary;

  public:
    Distributed(DBDataExchanger exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);

    template <typename DistributedBoundaryType>
    void ComputeBedFlux(DistributedBoundaryType& dbound);
};

Distributed::Distributed(DBDataExchanger exchanger) : exchanger(std::move(exchanger)) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {
    uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    this->q_ex.resize(SWE::n_variables, ngp);
    this->aux_ex.resize(SWE::n_auxiliaries, ngp);
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    const uint ngp = edge_dbound.edge_data.get_ngp();
    std::vector<double> message(1 + (SWE::n_variables + 1) * ngp);
    edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);

    const bool wet_in = edge_dbound.boundary.data.wet_dry_state.wet;
    const bool wet_ex = (bool)message[0];

    if (wet_in || wet_ex) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;
        auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex                            = ngp - gp - 1;
            this->aux_ex(SWE::Auxiliaries::bath, gp_ex) = message[1 + (SWE::n_variables + 1) * gp];
            for (uint var = 0; var < SWE::n_variables; ++var) {
                this->q_ex(var, gp_ex) = message[1 + (SWE::n_variables + 1) * gp + var + 1];
            }
        }
        row(this->aux_ex, SWE::Auxiliaries::h) =
            row(this->q_ex, SWE::Variables::ze) + row(this->aux_ex, SWE::Auxiliaries::bath);

        // Our definition of numerical flux implies q_hat = 0.5 * (q_in + q_ex)
        // Also assume bath_hat = 0.5*(bath_in +bath_ex)
        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath) =
            (row(boundary.aux_at_gp, SWE::Auxiliaries::bath) + row(this->aux_ex, SWE::Auxiliaries::bath)) / 2.0;
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
                    double temp_flux_hc = boundary.F_hat_at_gp(Variables::hc, gp);

                    SWE::get_Fn(column(edge_internal.q_hat_at_gp, gp),
                                column(edge_internal.aux_hat_at_gp, gp),
                                column(edge_dbound.boundary.surface_normal, gp),
                                column(boundary.F_hat_at_gp, gp),
                                0.0);

                    SWE::get_tau_LF(column(edge_internal.q_hat_at_gp, gp),
                                    column(edge_internal.aux_hat_at_gp, gp),
                                    column(edge_dbound.boundary.surface_normal, gp),
                                    edge_internal.tau[gp],
                                    0.0);

                    column(boundary.F_hat_at_gp, gp) +=
                        edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));

                    // Only remove gravity contributions for the momentum fluxes
                    boundary.F_hat_at_gp(Variables::ze, gp) = temp_flux_ze;
                    boundary.F_hat_at_gp(Variables::hc, gp) = temp_flux_hc;
                }
            }

            assert(!std::isnan(boundary.F_hat_at_gp(Variables::ze, gp)));
        }
    }
}

template <typename DistributedBoundaryType>
void Distributed::ComputeBedFlux(DistributedBoundaryType& dbound) {
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        const double un = roe_un(column(boundary.q_at_gp, gp),
                                 column(this->q_ex, gp),
                                 column(boundary.aux_at_gp, gp),
                                 column(this->aux_ex, gp),
                                 column(dbound.surface_normal, gp));
        if (Utilities::almost_equal(un, 0.0)) {
            boundary.qb_hat_at_gp[gp] = 0.0;
        } else if (un > 0.0) {
            boundary.qb_hat_at_gp[gp] = transpose(column(dbound.surface_normal, gp)) *
                                        bed_flux(column(boundary.q_at_gp, gp), column(boundary.aux_at_gp, gp));
        } else if (un < 0.0) {
            boundary.qb_hat_at_gp[gp] = transpose(column(dbound.surface_normal, gp)) *
                                        bed_flux(column(this->q_ex, gp), column(this->aux_ex, gp));
        }
    }
}
}
}
}

#endif