#ifndef SWE_SEABED_UPDATE_HPP
#define SWE_SEABED_UPDATE_HPP

#include "swe_seabed_CS_sl.hpp"
#include "swe_seabed_XU_sl.hpp"

namespace SWE {
StatVector<double, SWE::n_dimensions> bed_flux(const double h, const double qx, const double qy) {
    double usq = std::pow(qx / h, 2) + std::pow(qy / h, 2);
    double A   = 0.001;
    return StatVector<double, SWE::n_dimensions>{-A * usq * qx / h, -A * usq * qy / h};
}

double roe_un(const double nx,
              const double ny,
              const double h_in,
              const double qx_in,
              const double qy_in,
              const double h_ex,
              const double qx_ex,
              const double qy_ex) {
    return (qx_in / std::sqrt(h_in) + qx_ex / std::sqrt(h_ex)) / (std::sqrt(h_in) + std::sqrt(h_ex)) * nx +
           (qy_in / std::sqrt(h_in) + qy_ex / std::sqrt(h_ex)) / (std::sqrt(h_in) + std::sqrt(h_ex)) * ny;
}

template <typename StepperType, typename DiscretizationType>
void seabed_update(const StepperType& stepper, DiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& state    = elt.data.state[stepper.GetStage()];
        auto& internal = elt.data.internal;

        set_constant(state.b_rhs, 0.0);

        if (elt.data.wet_dry_state.wet) {
            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                column(internal.qb_at_gp, gp) =
                    bed_flux(internal.q_at_gp(SWE::Variables::ze, gp) + internal.aux_at_gp(SWE::Auxiliaries::bath, gp),
                             internal.q_at_gp(SWE::Variables::qx, gp),
                             internal.q_at_gp(SWE::Variables::qy, gp));
            }
            state.b_rhs = elt.IntegrationDPhi(GlobalCoord::x, row(internal.qb_at_gp, GlobalCoord::x)) +
                          elt.IntegrationDPhi(GlobalCoord::y, row(internal.qb_at_gp, GlobalCoord::y));
        }
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        if (intface.data_in.wet_dry_state.wet && intface.data_ex.wet_dry_state.wet) {
            intface.specialization.ComputeBedFlux(intface);
            state_in.b_rhs -= intface.IntegrationPhiIN(boundary_in.qb_hat_at_gp);
            state_ex.b_rhs -= intface.IntegrationPhiEX(boundary_ex.qb_hat_at_gp);
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& state    = bound.data.state[stepper.GetStage()];
        auto& boundary = bound.data.boundary[bound.bound_id];
        if (bound.data.wet_dry_state.wet) {
            bound.boundary_condition.ComputeBedFlux(stepper, bound);
            state.b_rhs -= bound.IntegrationPhi(boundary.qb_hat_at_gp);
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& state    = dbound.data.state[stepper.GetStage()];
        auto& boundary = dbound.data.boundary[dbound.bound_id];
        if (dbound.data.wet_dry_state.wet) {
            dbound.boundary_condition.ComputeBedFlux(dbound);
            state.b_rhs -= dbound.IntegrationPhi(boundary.qb_hat_at_gp);
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state;
        auto& next_state = elt.data.state[stage + 1];

        if (stage + 1 == stepper.GetNumStages()) {
            // swap back if we are at the last stage
            std::swap(state[0].q, state[stepper.GetNumStages()].q);
        }

        state[stage].b_solution = elt.ApplyMinv(state[stage].b_rhs);
        set_constant(row(next_state.aux, SWE::Auxiliaries::bath), 0.0);
        set_constant(row(next_state.q, SWE::Variables::ze), 0.0);
        for (uint s = 0; s <= stage; ++s) {
            row(next_state.aux, SWE::Auxiliaries::bath) +=
                stepper.ark[stage][s] * row(state[s].aux, SWE::Auxiliaries::bath) +
                stepper.GetDT() * stepper.brk[stage][s] * state[s].b_solution;

            row(next_state.q, SWE::Variables ::ze) +=
                stepper.ark[stage][s] * row(state[s].q, SWE::Variables::ze) +
                stepper.GetDT() * stepper.brk[stage][s] * row(state[s].solution, SWE::Variables::ze) -
                stepper.GetDT() * stepper.brk[stage][s] * state[s].b_solution;
        }

        if (stage + 1 == stepper.GetNumStages()) {
            std::swap(state[0].q, state[stepper.GetNumStages()].q);
            std::swap(state[0].aux, state[stepper.GetNumStages()].aux);
        }
    });
}

template <typename StepperType, typename DiscretizationType>
void seabed_data_update(const StepperType& stepper, DiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& state                                     = elt.data.state[stepper.GetStage()];
        auto& internal                                  = elt.data.internal;
        row(internal.aux_at_gp, SWE::Auxiliaries::bath) = elt.ComputeUgp(row(state.aux, SWE::Auxiliaries::bath));
        row(internal.db_at_gp, GlobalCoord::x) =
            elt.ComputeDUgp(GlobalCoord::x, row(state.aux, SWE::Auxiliaries::bath));
        row(internal.db_at_gp, GlobalCoord::y) =
            elt.ComputeDUgp(GlobalCoord::y, row(state.aux, SWE::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath) =
            intface.ComputeUgpIN(row(state_in.aux, SWE::Auxiliaries::bath));
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath) =
            intface.ComputeUgpEX(row(state_ex.aux, SWE::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& state                                     = bound.data.state[stepper.GetStage()];
        auto& boundary                                  = bound.data.boundary[bound.bound_id];
        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) = bound.ComputeUgp(row(state.aux, SWE::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& state                                     = dbound.data.state[stepper.GetStage()];
        auto& boundary                                  = dbound.data.boundary[dbound.bound_id];
        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) = dbound.ComputeUgp(row(state.aux, SWE::Auxiliaries::bath));
    });

    if (SWE::PostProcessing::wetting_drying) {
        discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            auto& state                   = elt.data.state[stepper.GetStage()];
            auto& wd_state                = elt.data.wet_dry_state;
            DynRowVector<double> bath_lin = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                wd_state.bath_at_vrtx[vrtx] = bath_lin[vrtx];
            }
            wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());
        });
    }

    if (SWE::PostProcessing::slope_limiting || SWE::PostProcessing::bed_slope_limiting) {
        discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            auto& state              = elt.data.state[stepper.GetStage()];
            auto& sl_state           = elt.data.slope_limit_state;
            sl_state.bath_lin        = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
            sl_state.bath_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.bath_lin);
            sl_state.bath_at_vrtx    = sl_state.bath_lin;
            sl_state.bath_at_midpts  = elt.ComputeLinearUmidpts(sl_state.bath_lin);
        });
    }
}
}

#endif