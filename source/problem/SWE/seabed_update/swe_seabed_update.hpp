#ifndef SWE_SEABED_UPDATE_HPP
#define SWE_SEABED_UPDATE_HPP

#include "swe_seabed_CS_sl.hpp"
#include "swe_seabed_XU_sl.hpp"

namespace SWE {
StatVector<double, SWE::n_dimensions> bed_flux(const double h, const double qx, const double qy) {
    double usq = std::pow(qx / h, 2) + std::pow(qy / h, 2);
    double A   = 0.005;
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
        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        set_constant(state.b_rhs, 0.0);

        if (elt.data.wet_dry_state.wet) {
            internal.q_at_gp = elt.ComputeUgp(state.q);
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
        auto& state_in    = intface.data_in.state[0];
        auto& state_ex    = intface.data_ex.state[0];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        if (intface.data_in.wet_dry_state.wet && intface.data_ex.wet_dry_state.wet) {
            boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
            boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);
            intface.specialization.ComputeBedFlux(intface);
            state_in.b_rhs -= intface.IntegrationPhiIN(boundary_in.qb_hat_at_gp);
            state_ex.b_rhs -= intface.IntegrationPhiEX(boundary_ex.qb_hat_at_gp);
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& state    = bound.data.state[0];
        auto& boundary = bound.data.boundary[bound.bound_id];
        if (bound.data.wet_dry_state.wet) {
            boundary.q_at_gp = bound.ComputeUgp(state.q);
            bound.boundary_condition.ComputeBedFlux(stepper, bound);
            state.b_rhs -= bound.IntegrationPhi(boundary.qb_hat_at_gp);
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& state    = dbound.data.state[0];
        auto& boundary = dbound.data.boundary[dbound.bound_id];
        if (dbound.data.wet_dry_state.wet) {
            boundary.q_at_gp = dbound.ComputeUgp(state.q);
            dbound.boundary_condition.ComputeBedFlux(dbound);
            state.b_rhs -= dbound.IntegrationPhi(boundary.qb_hat_at_gp);
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& state      = elt.data.state[0];
        state.b_solution = elt.ApplyMinv(state.b_rhs);
        row(state.aux, SWE::Auxiliaries::bath) += stepper.GetDT() * state.b_solution;
        row(state.q, SWE::Variables::ze) -= stepper.GetDT() * state.b_solution;
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

template <typename DiscretizationType>
void seabed_data_update(DiscretizationType& discretization) {
    // UPDATE VALUES
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;
        row(internal.aux_at_gp, SWE::Auxiliaries::bath) =
            elt.ComputeUgp(row(elt.data.state[0].aux, SWE::Auxiliaries::bath));
        row(internal.db_at_gp, GlobalCoord::x) =
            elt.ComputeDUgp(GlobalCoord::x, row(elt.data.state[0].aux, SWE::Auxiliaries::bath));
        row(internal.db_at_gp, GlobalCoord::y) =
            elt.ComputeDUgp(GlobalCoord::y, row(elt.data.state[0].aux, SWE::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath) =
            intface.ComputeUgpIN(row(intface.data_in.state[0].aux, SWE::Auxiliaries::bath));
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath) =
            intface.ComputeUgpEX(row(intface.data_ex.state[0].aux, SWE::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& boundary = bound.data.boundary[bound.bound_id];
        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) =
            bound.ComputeUgp(row(bound.data.state[0].aux, SWE::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& boundary = dbound.data.boundary[dbound.bound_id];
        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) =
            dbound.ComputeUgp(row(dbound.data.state[0].aux, SWE::Auxiliaries::bath));
    });

    if (SWE::PostProcessing::wetting_drying) {
        discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state                   = elt.data.state[0];
            auto& wd_state                = elt.data.wet_dry_state;
            DynRowVector<double> bath_lin = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                wd_state.bath_at_vrtx[vrtx] = bath_lin[vrtx];
            }
            wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());
        });
    }

    if (SWE::PostProcessing::slope_limiting || SWE::PostProcessing::bed_slope_limiting) {
        discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state              = elt.data.state[0];
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