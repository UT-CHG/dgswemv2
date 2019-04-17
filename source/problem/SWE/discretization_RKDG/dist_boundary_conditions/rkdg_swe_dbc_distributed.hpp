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
    std::array<DynRowVector<double>, SWE::n_variables> q_ex;

    BC::Land land_boundary;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    void ComputeFlux(DistributedBoundaryType& dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {
    uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    this->q_ex.fill(DynRowVector<double>(ngp));

    this->land_boundary.Initialize(1);
}

template <typename DistributedBoundaryType>
void Distributed::ComputeFlux(DistributedBoundaryType& dbound) {
    std::vector<double> message;

    message.resize(1 + SWE::n_variables * dbound.data.get_ngp_boundary(dbound.bound_id));

    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);

    bool wet_ex = (bool)message[0];

    uint gp_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        gp_ex = dbound.data.get_ngp_boundary(dbound.bound_id) - gp - 1;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            this->q_ex[var][gp_ex] = message[1 + SWE::n_variables * gp + var];
        }
    }

    bool wet_in = dbound.data.wet_dry_state.wet;

    auto& boundary = dbound.data.boundary[dbound.bound_id];

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        std::array<double, SWE::n_variables> F_hat_tmp;

        LLF_flux(Global::g,
                 boundary.q_at_gp[SWE::Variables::ze][gp],
                 boundary.q_at_gp[SWE::Variables::qx][gp],
                 boundary.q_at_gp[SWE::Variables::qy][gp],
                 this->q_ex[SWE::Variables::ze][gp],
                 this->q_ex[SWE::Variables::qx][gp],
                 this->q_ex[SWE::Variables::qy][gp],
                 std::array<double,SWE::n_auxiliaries>{boundary.aux_at_gp[SWE::Auxiliaries::bath][gp],
                         boundary.aux_at_gp[SWE::Auxiliaries::h][gp],
                         boundary.aux_at_gp[SWE::Auxiliaries::sp][gp]},
                 std::array<double,SWE::n_dimensions>{dbound.surface_normal[GlobalCoord::x][gp],
		     dbound.surface_normal[GlobalCoord::y][gp]},
                 F_hat_tmp
            );

	for ( uint i = 0; i < SWE::n_variables; ++i ) {
            boundary.F_hat_at_gp[i][gp] = F_hat_tmp[i];
        }

        /*if (boundary.F_hat_at_gp(Variables::ze, gp) > 1e-12) {
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

                LLF_flux(0.0,
                         column(boundary.q_at_gp, gp),
                         column(this->q_ex, gp),
                         column(boundary.aux_at_gp, gp),
                         column(dbound.surface_normal, gp),
                         column(boundary.F_hat_at_gp, gp));

                // Only remove gravity contributions for the momentum fluxes
                boundary.F_hat_at_gp(Variables::ze, gp) = temp_flux_ze;
            }
	    }

	    assert(!std::isnan(boundary.F_hat_at_gp(Variables::ze, gp)));*/
    }
}
}
}
}

#endif