#ifndef SWE_IS_LEVEE_HPP
#define SWE_IS_LEVEE_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/rkdg_simulation/rk_stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace IS {
class Levee {
  private:
    double H_tolerance = 0.01;

    std::vector<double> H_barrier;
    std::vector<double> C_subcritical;
    std::vector<double> C_supercritical;

    std::vector<double> H_bar_gp;
    std::vector<double> C_subcrit_gp;
    std::vector<double> C_supercrit_gp;

  public:
    Levee() = default;
    Levee(const std::vector<LeveeInput>& levee_input);

    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename InterfaceType>
    void ComputeFlux(const RKStepper& stepper, InterfaceType& intface);
};

Levee::Levee(const std::vector<LeveeInput>& levee_input) {
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

template <typename InterfaceType>
void Levee::Initialize(InterfaceType& intface) {
    this->H_bar_gp.resize(intface.data_in.get_ngp_boundary(intface.bound_id_in));
    this->C_subcrit_gp.resize(intface.data_in.get_ngp_boundary(intface.bound_id_in));
    this->C_supercrit_gp.resize(intface.data_in.get_ngp_boundary(intface.bound_id_in));

    intface.ComputeBoundaryNodalUgpIN(this->H_barrier, this->H_bar_gp);
    intface.ComputeBoundaryNodalUgpIN(this->C_subcritical, this->C_subcrit_gp);
    intface.ComputeBoundaryNodalUgpIN(this->C_supercritical, this->C_supercrit_gp);
}

template <typename InterfaceType>
void Levee::ComputeFlux(const RKStepper& stepper, InterfaceType& intface) {
    bool wet_in = intface.data_in.wet_dry_state.wet;
    bool wet_ex = intface.data_ex.wet_dry_state.wet;

    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& sp_at_gp_in = intface.data_in.spherical_projection.sp_at_gp_boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
    auto& sp_at_gp_ex = intface.data_ex.spherical_projection.sp_at_gp_boundary[intface.bound_id_ex];

    SWE::BC::Land land_boundary;

    double H_levee, C_subcrit, C_supercrit;
    double h_above_levee_in, h_above_levee_ex;
    double ze_in_ex, qx_in_ex, qy_in_ex, ze_ex_ex, qx_ex_ex, qy_ex_ex;
    double gravity_in, gravity_ex;

    uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);
    uint gp_ex = 0;
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        gp_ex = ngp - gp - 1;

        H_levee     = this->H_bar_gp[gp];
        C_subcrit   = this->C_subcrit_gp[gp];
        C_supercrit = this->C_supercrit_gp[gp];

        gravity_in = Global::g;
        gravity_ex = Global::g;

        h_above_levee_in = boundary_in.ze_at_gp[gp] - H_levee;
        h_above_levee_ex = boundary_ex.ze_at_gp[gp_ex] - H_levee;

        if ((h_above_levee_in <= H_tolerance && h_above_levee_ex <= H_tolerance) ||  // both side below or
            std::abs(h_above_levee_in - h_above_levee_ex) <= H_tolerance) {          // equal within tolerance
            // reflective boundary in
            land_boundary.GetEX(stepper,
                                gp,
                                intface.surface_normal_in,
                                boundary_in.ze_at_gp,
                                boundary_in.qx_at_gp,
                                boundary_in.qy_at_gp,
                                ze_in_ex,
                                qx_in_ex,
                                qy_in_ex);

            // reflective boundary ex
            land_boundary.GetEX(stepper,
                                gp_ex,
                                intface.surface_normal_ex,
                                boundary_ex.ze_at_gp,
                                boundary_ex.qx_at_gp,
                                boundary_ex.qy_at_gp,
                                ze_ex_ex,
                                qx_ex_ex,
                                qy_ex_ex);
        } else if (h_above_levee_in > h_above_levee_ex) {  // overtopping from in to ex
            double n_x, n_y, t_x, t_y, qn_in, qt_in;

            n_x = intface.surface_normal_in[gp][GlobalCoord::x];
            n_y = intface.surface_normal_in[gp][GlobalCoord::y];
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

            ze_in_ex = boundary_in.ze_at_gp[gp];
            qx_in_ex = qn_in * n_x + qt_in * t_x;
            qy_in_ex = qn_in * n_y + qt_in * t_y;

            ze_ex_ex = boundary_ex.ze_at_gp[gp_ex];
            qx_ex_ex = qx_in_ex;
            qy_ex_ex = qy_in_ex;

            if (!wet_ex) {
                gravity_ex = 0.0;
            }
        } else if (h_above_levee_in < h_above_levee_ex) {  // overtopping from ex to in
            double n_x, n_y, t_x, t_y, qn_ex, qt_ex;

            n_x = intface.surface_normal_ex[gp_ex][GlobalCoord::x];
            n_y = intface.surface_normal_ex[gp_ex][GlobalCoord::y];
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

            ze_ex_ex = boundary_ex.ze_at_gp[gp_ex];
            qx_ex_ex = qn_ex * n_x + qt_ex * t_x;
            qy_ex_ex = qn_ex * n_y + qt_ex * t_y;

            ze_in_ex = boundary_in.ze_at_gp[gp];
            qx_in_ex = qx_ex_ex;
            qy_in_ex = qy_ex_ex;

            if (!wet_in) {
                gravity_in = 0.0;
            }
        }

        LLF_flux(gravity_in,
                 boundary_in.ze_at_gp[gp],
                 ze_in_ex,
                 boundary_in.qx_at_gp[gp],
                 qx_in_ex,
                 boundary_in.qy_at_gp[gp],
                 qy_in_ex,
                 boundary_in.bath_at_gp[gp],
                 sp_at_gp_in[gp],
                 intface.surface_normal_in[gp],
                 boundary_in.ze_numerical_flux_at_gp[gp],
                 boundary_in.qx_numerical_flux_at_gp[gp],
                 boundary_in.qy_numerical_flux_at_gp[gp]);

        LLF_flux(gravity_ex,
                 boundary_ex.ze_at_gp[gp_ex],
                 ze_ex_ex,
                 boundary_ex.qx_at_gp[gp_ex],
                 qx_ex_ex,
                 boundary_ex.qy_at_gp[gp_ex],
                 qy_ex_ex,
                 boundary_ex.bath_at_gp[gp_ex],
                 sp_at_gp_ex[gp_ex],
                 intface.surface_normal_ex[gp_ex],
                 boundary_ex.ze_numerical_flux_at_gp[gp_ex],
                 boundary_ex.qx_numerical_flux_at_gp[gp_ex],
                 boundary_ex.qy_numerical_flux_at_gp[gp_ex]);
    }
}
}
}

#endif