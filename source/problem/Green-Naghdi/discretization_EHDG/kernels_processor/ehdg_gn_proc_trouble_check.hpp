#ifndef EHDG_GN_PROC_TROUBLE_CHECK_HPP
#define EHDG_GN_PROC_TROUBLE_CHECK_HPP

namespace GN {
namespace EHDG {
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

template <typename DiscretizationType, typename StepperType>
void check_trouble(DiscretizationType& discretization, const StepperType& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { elt.data.source.I = 0.0; });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];
        auto& source_in     = intface.data_in.source;
        auto& source_ex     = intface.data_ex.source;
        auto& derivative_in = intface.data_in.derivative;

        if ((intface.data_in.wet_dry_state.wet) && (intface.data_ex.wet_dry_state.wet)) {
            const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex = ngp - gp - 1;
                const double un  = intface.surface_normal_in(GlobalCoord::x, gp) *
                                      derivative_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::x, gp) +
                                  intface.surface_normal_in(GlobalCoord::y, gp) *
                                      derivative_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::y, gp);

                //                if (Utilities::almost_equal(un, 0.0)) {
                //                    boundary_in.hdif_at_gp[gp] = 0.0;
                //                    boundary_ex.hdif_at_gp[gp_ex] = 0.0;
                //                } else if (un < 0.0) {
                boundary_in.hdif_at_gp[gp] =
                    boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp) - boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex);
                //                    boundary_ex.hdif_at_gp[gp_ex] = 0.0;
                //                } else if (un > 0.0) {
                //                    boundary_in.hdif_at_gp[gp] = 0.0;
                boundary_ex.hdif_at_gp[gp_ex] =
                    boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex) - boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp);
                //                }
            }

            source_in.I += intface.IntegrationIN(boundary_in.hdif_at_gp);
            source_ex.I += intface.IntegrationEX(boundary_ex.hdif_at_gp);
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& boundary   = dbound.data.boundary[dbound.bound_id];
        auto& source     = dbound.data.source;
        auto& derivative = dbound.data.derivative;

        if (dbound.data.wet_dry_state.wet && dbound.boundary_condition.wet_neighbor) {
            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            for (uint gp = 0; gp < ngp; ++gp) {
                const double un = dbound.surface_normal(GlobalCoord::x, gp) *
                                      derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::x, gp) +
                                  dbound.surface_normal(GlobalCoord::y, gp) *
                                      derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::y, gp);

                //                if (Utilities::almost_equal(un, 0.0)) {
                //                   boundary.hdif_at_gp[gp] = 0.0;
                //                } else if (un < 0.0) {
                boundary.hdif_at_gp[gp] =
                    2.0 * (boundary.q_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp) -
                           boundary.aux_at_gp(SWE::Auxiliaries::h, gp));
                //                } else if (un > 0.0) {
                //                   boundary.hdif_at_gp[gp] = 0.0;
                //                }
            }

            source.I += dbound.Integration(boundary.hdif_at_gp);
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        if (elt.data.wet_dry_state.wet) {
            auto& state    = elt.data.state[stepper.GetStage()];
            auto& sl_state = elt.data.slope_limit_state;
            auto& source   = elt.data.source;

            sl_state.q_lin        = elt.ProjectBasisToLinear(state.q);
            sl_state.q_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.q_lin);
            const double h_norm   = sl_state.q_at_baryctr[SWE::Variables::ze] + sl_state.bath_at_baryctr;

            source.I /= (source.radius * source.perimeter * h_norm);
            if (elt.GetID() == 1440)
                std::cout << source.I << std::endl;
            if (std::abs(source.I) > 0.5) {
                elt.data.wet_dry_state.went_completely_dry = true;
            }
        }
    });
}
}
}

#endif