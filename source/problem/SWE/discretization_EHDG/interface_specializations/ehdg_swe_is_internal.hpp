#ifndef EHDG_SWE_IS_INTERNAL_HPP
#define EHDG_SWE_IS_INTERNAL_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace EHDG {
namespace ISP {
class Internal {
  private:
    BC::Land land_boundary;

  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFlux(EdgeInterfaceType& edge_int);

    template <typename InterfaceType>
    void ComputeBedFlux(InterfaceType& intface);
};

template <typename InterfaceType>
void Internal::Initialize(InterfaceType& intface) {}

template <typename EdgeInterfaceType>
void Internal::ComputeNumericalFlux(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;
    auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    const bool wet_in = edge_int.interface.data_in.wet_dry_state.wet;
    const bool wet_ex = edge_int.interface.data_ex.wet_dry_state.wet;

    // Our definition of numerical flux implies q_hat = 0.5 * (q_in + q_ex)
    const uint ngp = edge_int.edge_data.get_ngp();
    uint gp_ex;
    for (uint gp = 0; gp < ngp; ++gp) {
        gp_ex = ngp - gp - 1;

        column(edge_internal.q_hat_at_gp, gp) =
            (column(boundary_in.q_at_gp, gp) + column(boundary_ex.q_at_gp, gp_ex)) / 2.0;
    }

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    /* Compute trace flux for in side */

    SWE::get_Fn(edge_internal.q_hat_at_gp,
                edge_internal.aux_hat_at_gp,
                edge_int.interface.surface_normal_in,
                boundary_in.F_hat_at_gp);

    /* Add stabilization parameter terms */

    SWE::get_tau_LF(edge_internal.q_hat_at_gp,
                    edge_internal.aux_hat_at_gp,
                    edge_int.interface.surface_normal_in,
                    edge_internal.tau);

    gp_ex = 0;
    for (uint gp = 0; gp < ngp; ++gp) {
        gp_ex = ngp - gp - 1;

        column(boundary_ex.F_hat_at_gp, gp) = -column(boundary_in.F_hat_at_gp, gp_ex);

        column(boundary_in.F_hat_at_gp, gp) +=
            edge_internal.tau[gp] * (column(boundary_in.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
        column(boundary_ex.F_hat_at_gp, gp_ex) +=
            edge_internal.tau[gp] * (column(boundary_ex.q_at_gp, gp_ex) - column(edge_internal.q_hat_at_gp, gp));
    }

    // corrent numerical fluxes for wetting/drying
    gp_ex = 0;
    for (uint gp = 0; gp < ngp; ++gp) {
        gp_ex = ngp - gp - 1;

        if (boundary_in.F_hat_at_gp(Variables::ze, gp) > 1e-12) {
            if (!wet_in) {  // water flowing from dry IN element
                // Zero flux on IN element side
                set_constant(column(boundary_in.F_hat_at_gp, gp), 0.0);

                // Reflective Boundary on EX element side
                this->land_boundary.ComputeFlux(column(edge_int.interface.surface_normal_ex, gp_ex),
                                                column(boundary_ex.q_at_gp, gp_ex),
                                                column(edge_internal.aux_hat_at_gp, gp),
                                                edge_internal.tau[gp],
                                                column(edge_internal.q_hat_at_gp, gp),
                                                column(boundary_ex.F_hat_at_gp, gp_ex));

            } else if (!wet_ex) {  // water flowing to dry EX element
                SWE::get_Fn(0.0,
                            column(edge_internal.q_hat_at_gp, gp),
                            column(edge_internal.aux_hat_at_gp, gp),
                            column(edge_int.interface.surface_normal_ex, gp_ex),
                            column(boundary_ex.F_hat_at_gp, gp_ex));

                SWE::get_tau_LF(0.0,
                                column(edge_internal.q_hat_at_gp, gp),
                                column(edge_internal.aux_hat_at_gp, gp),
                                column(edge_int.interface.surface_normal_ex, gp_ex),
                                edge_internal.tau[gp]);

                column(boundary_ex.F_hat_at_gp, gp_ex) +=
                    edge_internal.tau[gp] *
                    (column(boundary_ex.q_at_gp, gp_ex) - column(edge_internal.q_hat_at_gp, gp));

                // Only remove gravity contributions for the momentum fluxes
                boundary_ex.F_hat_at_gp(Variables::ze, gp_ex) = -boundary_in.F_hat_at_gp(Variables::ze, gp);
            }
        } else if (boundary_in.F_hat_at_gp(Variables::ze, gp) < -1e-12) {
            if (!wet_ex) {  // water flowing from dry EX element
                // Zero flux on EX element side
                set_constant(column(boundary_ex.F_hat_at_gp, gp_ex), 0.0);

                // Reflective Boundary on IN element side

                this->land_boundary.ComputeFlux(column(edge_int.interface.surface_normal_in, gp),
                                                column(boundary_in.q_at_gp, gp),
                                                column(edge_internal.aux_hat_at_gp, gp),
                                                edge_internal.tau[gp],
                                                column(edge_internal.q_hat_at_gp, gp),
                                                column(boundary_in.F_hat_at_gp, gp));

            } else if (!wet_in) {  // water flowing to dry IN element
                SWE::get_Fn(0.0,
                            column(edge_internal.q_hat_at_gp, gp),
                            column(edge_internal.aux_hat_at_gp, gp),
                            column(edge_int.interface.surface_normal_in, gp),
                            column(boundary_in.F_hat_at_gp, gp));

                SWE::get_tau_LF(0.0,
                                column(edge_internal.q_hat_at_gp, gp),
                                column(edge_internal.aux_hat_at_gp, gp),
                                column(edge_int.interface.surface_normal_in, gp),
                                edge_internal.tau[gp]);

                column(boundary_in.F_hat_at_gp, gp) +=
                    edge_internal.tau[gp] * (column(boundary_in.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));

                // Only remove gravity contributions for the momentum fluxes
                boundary_in.F_hat_at_gp(Variables::ze, gp) = -boundary_ex.F_hat_at_gp(Variables::ze, gp_ex);
            }
        }

        assert(!std::isnan(boundary_in.F_hat_at_gp(Variables::ze, gp)));
    }
}

template <typename InterfaceType>
void Internal::ComputeBedFlux(InterfaceType& intface) {
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
    for (uint gp = 0; gp < ngp; ++gp) {
        const uint gp_ex = ngp - gp - 1;
        const double un  = roe_un(column(boundary_in.q_at_gp, gp),
                                 column(boundary_ex.q_at_gp, gp_ex),
                                 column(boundary_in.aux_at_gp, gp),
                                 column(boundary_ex.aux_at_gp, gp_ex),
                                 column(intface.surface_normal_in, gp));

        if (Utilities::almost_equal(un, 0.0)) {
            boundary_in.qb_hat_at_gp[gp]    = 0.0;
            boundary_ex.qb_hat_at_gp[gp_ex] = 0.0;
        } else if (un > 0.0) {
            boundary_in.qb_hat_at_gp[gp] = transpose(column(intface.surface_normal_in, gp)) *
                                           bed_flux(column(boundary_in.q_at_gp, gp), column(boundary_in.aux_at_gp, gp));
            boundary_ex.qb_hat_at_gp[gp_ex] = -boundary_in.qb_hat_at_gp[gp];
        } else if (un < 0.0) {
            boundary_in.qb_hat_at_gp[gp] =
                transpose(column(intface.surface_normal_in, gp)) *
                bed_flux(column(boundary_ex.q_at_gp, gp_ex), column(boundary_ex.aux_at_gp, gp_ex));
            boundary_ex.qb_hat_at_gp[gp_ex] = -boundary_in.qb_hat_at_gp[gp];
        }
    }
}
}
}
}

#endif