#ifndef RKDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP
#define RKDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

#include "rkdg_swe_dbc_distributed.hpp"

namespace SWE {
namespace RKDG {
namespace DBC {
class DistributedLevee {
  public:
    DBDataExchanger exchanger;

    double H_tolerance = 0.01;

    std::vector<double> H_barrier;
    std::vector<double> C_subcritical;
    std::vector<double> C_supercritical;

    std::vector<double> H_bar_gp;
    std::vector<double> C_subcrit_gp;
    std::vector<double> C_supercrit_gp;

  public:
    DistributedLevee() = default;
    DistributedLevee(const DBDataExchanger& exchanger, const std::vector<LeveeInput>& levee_input);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    void ComputeFlux(const RKStepper& stepper, DistributedBoundaryType& dbound);
};

DistributedLevee::DistributedLevee(const DBDataExchanger& exchanger, const std::vector<LeveeInput>& levee_input)
    : exchanger(exchanger) {
    uint n_nodes = levee_input.size();

    H_barrier.resize(n_nodes);
    C_subcritical.resize(n_nodes);
    C_supercritical.resize(n_nodes);

    for (uint node = 0; node < n_nodes; node++) {
        this->H_barrier[node]       = levee_input[node].H_barrier;
        this->C_subcritical[node]   = levee_input[node].C_subcritical;
        this->C_supercritical[node] = levee_input[node].C_supercritical;
    }
}

template <typename DistributedBoundaryType>
void DistributedLevee::Initialize(DistributedBoundaryType& dbound) {
    this->H_bar_gp.resize(dbound.data.get_ngp_boundary(dbound.bound_id));
    this->C_subcrit_gp.resize(dbound.data.get_ngp_boundary(dbound.bound_id));
    this->C_supercrit_gp.resize(dbound.data.get_ngp_boundary(dbound.bound_id));

    dbound.ComputeBoundaryNodalUgp(this->H_barrier, this->H_bar_gp);
    dbound.ComputeBoundaryNodalUgp(this->C_subcritical, this->C_subcrit_gp);
    dbound.ComputeBoundaryNodalUgp(this->C_supercritical, this->C_supercrit_gp);
}

template <typename DistributedBoundaryType>
void DistributedLevee::ComputeFlux(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    bool wet_in = dbound.data.wet_dry_state.wet;
    bool wet_ex;
    this->exchanger.GetWetDryEX(wet_ex);

    auto& boundary = dbound.data.boundary[dbound.bound_id];

    BC::Land land_boundary;

    double H_levee, C_subcrit, C_supercrit;
    double h_above_levee_in, h_above_levee_ex;
    StatVector<double, SWE::n_variables> q_ex;
    double gravity;

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        this->exchanger.GetEX(gp, q_ex);

        H_levee     = this->H_bar_gp[gp];
        C_subcrit   = this->C_subcrit_gp[gp];
        C_supercrit = this->C_supercrit_gp[gp];

        gravity = Global::g;

        h_above_levee_in = boundary.q_at_gp(SWE::Variables::ze, gp) - H_levee;
        h_above_levee_ex = q_ex[SWE::Variables::ze] - H_levee;

        if ((h_above_levee_in <= H_tolerance && h_above_levee_ex <= H_tolerance) ||  // both side below or
            std::abs(h_above_levee_in - h_above_levee_ex) <= H_tolerance) {          // equal within tolerance
            // reflective boundary
            land_boundary.GetEX(stepper, dbound.surface_normal[gp], boundary.q_at_gp[gp], q_ex);
        } else if (h_above_levee_in > h_above_levee_ex) {  // overtopping from in to ex
            double n_x, n_y, t_x, t_y, qn_in, qt_in;

            n_x = dbound.surface_normal[gp][GlobalCoord::x];
            n_y = dbound.surface_normal[gp][GlobalCoord::y];
            t_x = -n_y;
            t_y = n_x;

            if (h_above_levee_ex > 2.0 * h_above_levee_in / 3.0) {  // subcritical flow
                qn_in =
                    C_subcrit * h_above_levee_ex * std::sqrt(2.0 * Global::g * (h_above_levee_in - h_above_levee_ex));
                qt_in = 0.0;
            } else {  // supercritical flow
                qn_in = C_supercrit * (2.0 * h_above_levee_in / 3.0) *
                        std::sqrt(Global::g * (2.0 * h_above_levee_in / 3.0));
                qt_in = 0.0;
            }

            q_ex[SWE::Variables::ze] = boundary.q_at_gp(SWE::Variables::ze, gp);
            q_ex[SWE::Variables::qx] = qn_in * n_x + qt_in * t_x;
            q_ex[SWE::Variables::qy] = qn_in * n_y + qt_in * t_y;
        } else if (h_above_levee_in < h_above_levee_ex) {  // overtopping from ex to in
            double n_x, n_y, t_x, t_y, qn_ex, qt_ex;

            n_x = dbound.surface_normal[gp][GlobalCoord::x];
            n_y = dbound.surface_normal[gp][GlobalCoord::y];
            t_x = -n_y;
            t_y = n_x;

            if (h_above_levee_in > 2.0 * h_above_levee_ex / 3.0) {  // subcritical flow
                qn_ex =
                    C_subcrit * h_above_levee_in * std::sqrt(2.0 * Global::g * (h_above_levee_ex - h_above_levee_in));
                qt_ex = 0.0;
            } else {  // supercritical flow
                qn_ex = C_supercrit * (2.0 * h_above_levee_ex / 3.0) *
                        std::sqrt(Global::g * (2.0 * h_above_levee_ex / 3.0));
                qt_ex = 0.0;
            }

            q_ex[SWE::Variables::ze] = boundary.q_at_gp(SWE::Variables::ze, gp);
            q_ex[SWE::Variables::qx] = -(qn_ex * n_x + qt_ex * t_x);
            q_ex[SWE::Variables::qy] = -(qn_ex * n_y + qt_ex * t_y);

            if (!wet_in) {
                gravity = 0.0;
            }
        }

        column(boundary.F_hat_at_gp, gp) = LLF_flux(
            gravity, column(boundary.q_at_gp, gp), q_ex, column(boundary.aux_at_gp, gp), dbound.surface_normal[gp]);
    }
}
}
}
}

#endif