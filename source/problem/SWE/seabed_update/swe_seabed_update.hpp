#ifndef SWE_SEABED_UPDATE_HPP
#define SWE_SEABED_UPDATE_HPP

#include "swe_seabed_CS_sl.hpp"
#include "swe_seabed_XU_sl.hpp"

namespace SWE {
double rho_mixture(const Column<HybMatrix<double, SWE::n_variables>>& q,
                   const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux) {
    const double c = q[SWE::Variables::hc] / aux[SWE::Auxiliaries::h];
    return SWE::Global::rho_water * (1.0 - c) + SWE::Global::rho_sediment * c;
}

double entrainment_rate(const Column<HybMatrix<double, SWE::n_variables>>& q,
                        const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux,
                        const bool manning,
                        const double gn_sq) {
    if (SWE::SedimentTransport::suspended_load) {
        constexpr double phi     = 0.01;
        constexpr double d       = 4.0e-3;
        constexpr double theta_c = 0.045;

        const double h = aux[SWE::Auxiliaries::h];
        const double u = std::hypot(q[SWE::Variables::qx], q[SWE::Variables::qy]) / aux[SWE::Auxiliaries::h];
        const double s = SWE::Global::rho_sediment / SWE::Global::rho_water - 1.0;

        double Cf = SWE::SourceTerms::Cf;
        if (manning) {
            Cf = gn_sq / std::pow(h, 1.0 / 3.0);
            if (Cf < SWE::SourceTerms::Cf)
                Cf = SWE::SourceTerms::Cf;
        }
        const double theta = Cf * std::pow(u, 2) / (s * SWE::Global::g * d);

        if (theta > theta_c) {
            return phi * (theta - theta_c) * std::pow(d, -0.2) * u / h;
        }
    }
    return 0.0;
}

double deposition_rate(const Column<HybMatrix<double, SWE::n_variables>>& q,
                       const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux) {
    if (SWE::SedimentTransport::suspended_load) {
        constexpr double d  = 4.0e-3;
        constexpr double nu = 1.2e-6;

        const double c     = q[SWE::Variables::hc] / aux[SWE::Auxiliaries::h];
        const double s     = SWE::Global::rho_sediment / SWE::Global::rho_water - 1.0;
        const double alpha = std::min(2.0, (1.0 - SWE::Global::sat_sediment) / c);
        const double w_o   = std::sqrt(std::pow(13.95 * nu / d, 2) + 1.09 * s * Global::g * d) - 13.95 * nu / d;

        if (c > 0) {
            return alpha * c * w_o * std::pow(1 - alpha * c, 2);
        }
    }
    return 0.0;
}

StatVector<double, SWE::n_dimensions> bed_flux(const Column<HybMatrix<double, SWE::n_variables>>& q,
                                               const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux) {
    if (SWE::SedimentTransport::bed_load) {
        /*const double theta_c_o = 0.05;
        const double mu_s      = 0.63;
        const double mu_d      = 0.51;
        const double s         = 2.6;
        const double d         = 0.00018;

        const double alpha = std::acos((qx * bx + qy * by) / (std::hypot(qx, qy) * std::hypot(bx, by)));
        const double beta  = std::atan(std::hypot(bx, by));
        const double theta_c =
            theta_c_o * (std::cos(beta) *
                            std::sqrt(1 - std::pow(std::sin(alpha), 2) * std::pow(std::tan(beta), 2) / std::pow(mu_s,
        2)) - std::cos(alpha) * std::sin(beta) / mu_s); const double theta = std::pow(std::hypot(qx, qy) / h, 2) / ((s -
        1) * Global::g * d); const double p     = std::pow(1 + std::pow(PI * mu_d / (6 * (theta - theta_c)), 4), -0.25);
        const double A     = 100 * PI * d * p / 6;*/

        const double h  = aux[SWE::Auxiliaries::h];
        const double qx = q[SWE::Variables::qx];
        const double qy = q[SWE::Variables::qy];

        const double usq = std::pow(qx / h, 2) + std::pow(qy / h, 2);
        const double A   = 0.005;

        return StatVector<double, SWE::n_dimensions>{-A * usq * qx / h, -A * usq * qy / h};
    }
    return StatVector<double, SWE::n_dimensions>{0.0, 0.0};
}

double roe_un(const Column<HybMatrix<double, SWE::n_variables>>& q_in,
              const Column<HybMatrix<double, SWE::n_variables>>& q_ex,
              const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
              const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_ex,
              const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal) {
    const double h_in  = aux_in[SWE::Auxiliaries::h];
    const double qx_in = q_in[SWE::Variables::qx];
    const double qy_in = q_in[SWE::Variables::qy];

    const double h_ex  = aux_ex[SWE::Auxiliaries::h];
    const double qx_ex = q_ex[SWE::Variables::qx];
    const double qy_ex = q_ex[SWE::Variables::qy];

    const double nx = surface_normal[GlobalCoord::x];
    const double ny = surface_normal[GlobalCoord::y];

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
                column(internal.qb_at_gp, gp) = bed_flux(column(internal.q_at_gp, gp), column(internal.aux_at_gp, gp));
            }
            state.b_rhs = elt.IntegrationDPhi(GlobalCoord::x, row(internal.qb_at_gp, GlobalCoord::x)) +
                          elt.IntegrationDPhi(GlobalCoord::y, row(internal.qb_at_gp, GlobalCoord::y));

            state.b_rhs +=
                elt.IntegrationPhi(1.0 / (1.0 - Global::sat_sediment) * (internal.E_at_gp - internal.D_at_gp));
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

    /*std::set<uint> nodeIDs;
    discretization.mesh.CallForEachElement(
        [&nodeIDs](auto& elt) { nodeIDs.insert(elt.GetNodeID().begin(), elt.GetNodeID().end()); });
    uint max_nodeID = *std::max_element(nodeIDs.begin(), nodeIDs.end());

    DynVector<double> bath_at_node(max_nodeID + 1);
    set_constant(bath_at_node, 0.0);
    std::vector<uint> node_mult((max_nodeID + 1), 0);

    discretization.mesh.CallForEachElement([&bath_at_node, &node_mult](auto& elt) {
        auto& state    = elt.data.state[0];
        auto& sl_state = elt.data.slope_limit_state;

        sl_state.bath_lin = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
        sl_state.q_lin    = elt.ProjectBasisToLinear(state.q);
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            bath_at_node[elt.GetNodeID()[node]] += sl_state.bath_lin[node];
            ++node_mult[elt.GetNodeID()[node]];
        }
    });

    discretization.mesh.CallForEachElement([&bath_at_node, &node_mult](auto& elt) {
        auto& state    = elt.data.state[0];
        auto& sl_state = elt.data.slope_limit_state;

        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            sl_state.bath_at_vrtx[node] = bath_at_node[elt.GetNodeID()[node]] / node_mult[elt.GetNodeID()[node]];
        }
        row(elt.data.state[0].aux, SWE::Auxiliaries::bath) = elt.ProjectLinearToBasis(sl_state.bath_at_vrtx);
        row(sl_state.q_lin, SWE::Variables::ze) += (sl_state.bath_lin - sl_state.bath_at_vrtx);
        row(elt.data.state[0].q, SWE::Variables::ze) =
            elt.ProjectLinearToBasis(row(sl_state.q_lin, SWE::Variables::ze));
    });*/
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

    if (SWE::PostProcessing::slope_limiting || SWE::SedimentTransport::bed_slope_limiting) {
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