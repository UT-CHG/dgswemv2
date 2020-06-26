#ifndef RKDG_SWE_DBC_DISTRIBUTED_HPP
#define RKDG_SWE_DBC_DISTRIBUTED_HPP

#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace RKDG {
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

    template <typename DistributedBoundaryType>
    void ComputeFlux(DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    void ComputeBedFlux(DistributedBoundaryType& dbound);
};

Distributed::Distributed(DBDataExchanger exchanger) : exchanger(std::move(exchanger)) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {
    const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    this->q_ex.resize(SWE::n_variables, ngp);
    this->aux_ex.resize(SWE::n_auxiliaries, ngp);
    this->land_boundary.Initialize(1);
}

template <typename DistributedBoundaryType>
void Distributed::ComputeFlux(DistributedBoundaryType& dbound) {
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    std::vector<double> message(1 + (SWE::n_variables + 1) * ngp);
    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);
    for (uint gp = 0; gp < ngp; ++gp) {
        const uint gp_ex                            = ngp - gp - 1;
        this->aux_ex(SWE::Auxiliaries::bath, gp_ex) = message[1 + (SWE::n_variables + 1) * gp];
        for (uint var = 0; var < SWE::n_variables; ++var) {
            this->q_ex(var, gp_ex) = message[1 + (SWE::n_variables + 1) * gp + var + 1];
        }
    }
    row(this->aux_ex, SWE::Auxiliaries::h) =
        row(this->q_ex, SWE::Variables::ze) + row(this->aux_ex, SWE::Auxiliaries::bath);

    const bool wet_in = dbound.data.wet_dry_state.wet;
    const bool wet_ex = (bool)message[0];

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        HLL_flux(Global::g,
                 column(boundary.q_at_gp, gp),
                 column(this->q_ex, gp),
                 column(boundary.aux_at_gp, gp),
                 column(this->aux_ex, gp),
                 column(dbound.surface_normal, gp),
                 column(boundary.F_hat_at_gp, gp));

        if (boundary.F_hat_at_gp(Variables::ze, gp) > 1e-12) {
            if (!wet_in) {  // water flowing from dry IN element
                // Zero flux on IN element side
                set_constant(column(boundary.F_hat_at_gp, gp), 0.0);
            }
        } else if (boundary.F_hat_at_gp(Variables::ze, gp) < -1e-12) {
            if (!wet_ex) {  // water flowing from dry EX element
                // Reflective Boundary on IN element side
                this->land_boundary.ComputeFlux(column(dbound.surface_normal, gp),
                                                column(boundary.q_at_gp, gp),
                                                column(boundary.aux_at_gp, gp),
                                                column(boundary.F_hat_at_gp, gp));
            } else if (!wet_in) {  // water flowing to dry IN element
                double temp_flux_ze = boundary.F_hat_at_gp(Variables::ze, gp);
                double temp_flux_hc = boundary.F_hat_at_gp(Variables::hc, gp);

                HLL_flux(0.0,
                         column(boundary.q_at_gp, gp),
                         column(this->q_ex, gp),
                         column(boundary.aux_at_gp, gp),
                         column(this->aux_ex, gp),
                         column(dbound.surface_normal, gp),
                         column(boundary.F_hat_at_gp, gp));

                // Only remove gravity contributions for the momentum fluxes
                boundary.F_hat_at_gp(Variables::ze, gp) = temp_flux_ze;
                boundary.F_hat_at_gp(Variables::hc, gp) = temp_flux_hc;
            }
        }

        assert(!std::isnan(boundary.F_hat_at_gp(Variables::ze, gp)));
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