#ifndef SWE_INTFACE_SPEC_LEVEE_HPP
#define SWE_INTFACE_SPEC_LEVEE_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace IS {
class Levee {
  public:
    template <typename InterfaceType>
    void ComputeFlux(const Stepper& stepper, InterfaceType& intface);
};

template <typename InterfaceType>
void Levee::ComputeFlux(const Stepper& stepper, InterfaceType& intface) {
    double H_levee     = 0.35;
    double H_tolerance = 0.01;

    double C_subcrit   = 1.0;
    double C_supercrit = 1.0;

    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& sp_at_gp_in = intface.data_in.spherical_projection.sp_at_gp_boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
    auto& sp_at_gp_ex = intface.data_ex.spherical_projection.sp_at_gp_boundary[intface.bound_id_ex];

    SWE::BC::Land land_boundary;

    double h_above_levee_in;
    double h_above_levee_ex;
    double ze_in_ex, qx_in_ex, qy_in_ex, ze_ex_ex, qx_ex_ex, qy_ex_ex;

    uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);
    uint gp_ex = 0;
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        gp_ex = ngp - gp - 1;

        h_above_levee_in = boundary_in.ze_at_gp[gp] - H_levee;
        h_above_levee_ex = boundary_ex.ze_at_gp[gp_ex] - H_levee;

        if ((h_above_levee_in <= H_tolerance && h_above_levee_ex <= H_tolerance) ||
            std::abs(h_above_levee_in - h_above_levee_ex) <= H_tolerance) {  // both side below or equal within
                                                                             // tolerance
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
        }

        LLF_flux(Global::g,
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

        LLF_flux(Global::g,
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