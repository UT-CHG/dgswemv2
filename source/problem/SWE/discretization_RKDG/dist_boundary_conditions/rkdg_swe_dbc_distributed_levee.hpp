#ifndef RKDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP
#define RKDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP

#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace RKDG {
namespace DBC {
class DistributedLevee {
  public:
    DBDataExchanger exchanger;

  private:
    double H_tolerance = 0.01;

    HybMatrix<double, SWE::n_variables> q_ex;

    DynRowVector<double> H_barrier;
    DynRowVector<double> C_subcritical;
    DynRowVector<double> C_supercritical;

    DynRowVector<double> H_bar_gp;
    DynRowVector<double> C_subcrit_gp;
    DynRowVector<double> C_supercrit_gp;

    BC::Land land_boundary;

  public:
    DistributedLevee(DBDataExchanger  exchanger, const std::vector<LeveeInput>& levee_input);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    void ComputeFlux(DistributedBoundaryType& dbound);
};

DistributedLevee::DistributedLevee(DBDataExchanger  exchanger, const std::vector<LeveeInput>& levee_input)
    : exchanger(std::move(exchanger)) {
    uint n_nodes = levee_input.size();

    this->H_barrier.resize(n_nodes);
    this->C_subcritical.resize(n_nodes);
    this->C_supercritical.resize(n_nodes);

    for (uint node = 0; node < n_nodes; ++node) {
        this->H_barrier[node]       = levee_input[node].H_barrier;
        this->C_subcritical[node]   = levee_input[node].C_subcritical;
        this->C_supercritical[node] = levee_input[node].C_supercritical;
    }
}

template <typename DistributedBoundaryType>
void DistributedLevee::Initialize(DistributedBoundaryType& dbound) {
    uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

    this->q_ex.resize(SWE::n_variables, ngp);

    this->H_bar_gp       = dbound.ComputeBoundaryNodalUgp(this->H_barrier);
    this->C_subcrit_gp   = dbound.ComputeBoundaryNodalUgp(this->C_subcritical);
    this->C_supercrit_gp = dbound.ComputeBoundaryNodalUgp(this->C_supercritical);
}

template <typename DistributedBoundaryType>
void DistributedLevee::ComputeFlux(DistributedBoundaryType& dbound) {
    std::vector<double> message;

    message.resize(1 + SWE::n_variables * dbound.data.get_ngp_boundary(dbound.bound_id));

    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);

    uint gp_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        gp_ex = dbound.data.get_ngp_boundary(dbound.bound_id) - gp - 1;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            this->q_ex(var, gp_ex) = message[1 + SWE::n_variables * gp + var];
        }
    }

    bool wet_in = dbound.data.wet_dry_state.wet;

    auto& boundary = dbound.data.boundary[dbound.bound_id];

    double H_levee, C_subcrit, C_supercrit;
    double h_above_levee_in, h_above_levee_ex;
    double gravity;

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        H_levee     = this->H_bar_gp[gp];
        C_subcrit   = this->C_subcrit_gp[gp];
        C_supercrit = this->C_supercrit_gp[gp];

        gravity = Global::g;

        h_above_levee_in = boundary.q_at_gp(SWE::Variables::ze, gp) - H_levee;
        h_above_levee_ex = this->q_ex(SWE::Variables::ze, gp) - H_levee;

        if ((h_above_levee_in <= H_tolerance && h_above_levee_ex <= H_tolerance) ||  // both side below or
            std::abs(h_above_levee_in - h_above_levee_ex) <= H_tolerance) {          // equal within tolerance

            // reflective boundary
            this->land_boundary.GetEX(
                column(dbound.surface_normal, gp), column(boundary.q_at_gp, gp), column(this->q_ex, gp));
        } else if (h_above_levee_in > h_above_levee_ex) {  // overtopping from in to ex
            double n_x, n_y, t_x, t_y, qn_in, qt_in;

            n_x = dbound.surface_normal(GlobalCoord::x, gp);
            n_y = dbound.surface_normal(GlobalCoord::y, gp);
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

            this->q_ex(SWE::Variables::ze, gp) = boundary.q_at_gp(SWE::Variables::ze, gp);
            this->q_ex(SWE::Variables::qx, gp) = qn_in * n_x + qt_in * t_x;
            this->q_ex(SWE::Variables::qy, gp) = qn_in * n_y + qt_in * t_y;
        } else if (h_above_levee_in < h_above_levee_ex) {  // overtopping from ex to in
            double n_x, n_y, t_x, t_y, qn_ex, qt_ex;

            n_x = dbound.surface_normal(GlobalCoord::x, gp);
            n_y = dbound.surface_normal(GlobalCoord::y, gp);
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

            this->q_ex(SWE::Variables::ze, gp) = boundary.q_at_gp(SWE::Variables::ze, gp);
            this->q_ex(SWE::Variables::qx, gp) = -(qn_ex * n_x + qt_ex * t_x);
            this->q_ex(SWE::Variables::qy, gp) = -(qn_ex * n_y + qt_ex * t_y);

            if (!wet_in) {
                gravity = 0.0;
            }
        }

        LLF_flux(gravity,
                 column(boundary.q_at_gp, gp),
                 column(q_ex, gp),
                 column(boundary.aux_at_gp, gp),
                 column(dbound.surface_normal, gp),
                 column(boundary.F_hat_at_gp, gp));
    }
}
}
}
}

#endif