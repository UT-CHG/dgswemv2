#ifndef SWE_SEABED_CS_SL_HPP
#define SWE_SEABED_CS_SL_HPP

namespace SWE {
template <typename DiscretizationType>
void CS_seabed_slope_limiter(DiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state    = elt.data.state[0];
        auto& wd_state = elt.data.wet_dry_state;
        auto& sl_state = elt.data.slope_limit_state;
        if (wd_state.wet) {
            sl_state.bath_lin        = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
            sl_state.bath_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.bath_lin);
            sl_state.bath_at_midpts  = elt.ComputeLinearUmidpts(sl_state.bath_lin);
            sl_state.q_lin           = elt.ProjectBasisToLinear(state.q);
            sl_state.q_at_baryctr    = elt.ComputeLinearUbaryctr(sl_state.q_lin);
            sl_state.q_at_midpts     = elt.ComputeLinearUmidpts(sl_state.q_lin);
        }
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& wd_state_in = intface.data_in.wet_dry_state;
        auto& wd_state_ex = intface.data_ex.wet_dry_state;
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.wet_neigh[intface.bound_id_in] = wd_state_ex.wet;
        sl_state_ex.wet_neigh[intface.bound_id_ex] = wd_state_in.wet;
        if (wd_state_in.wet && wd_state_ex.wet) {
            sl_state_in.bath_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.bath_at_baryctr;
            sl_state_ex.bath_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.bath_at_baryctr;
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& wd_state = bound.data.wet_dry_state;
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.wet_neigh[bound.bound_id] = wd_state.wet;
        if (wd_state.wet) {
            sl_state.bath_at_baryctr_neigh[bound.bound_id] = sl_state.bath_at_baryctr;
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& wd_state = elt.data.wet_dry_state;
        auto& sl_state = elt.data.slope_limit_state;

        if (wd_state.wet &&
            std::find(sl_state.wet_neigh.begin(), sl_state.wet_neigh.end(), false) == sl_state.wet_neigh.end()) {
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                const uint element_1 = sl_state.a_elem[2 * bound];
                const uint element_2 = sl_state.a_elem[2 * bound + 1];

                const double b_tilda = sl_state.bath_at_midpts[bound] - sl_state.bath_at_baryctr;
                const double b_delta =
                    sl_state.alpha[bound][0] * (sl_state.bath_at_baryctr_neigh[element_1] - sl_state.bath_at_baryctr) +
                    sl_state.alpha[bound][1] * (sl_state.bath_at_baryctr_neigh[element_2] - sl_state.bath_at_baryctr);
                if (std::abs(b_tilda) <= PostProcessing::M * sl_state.r_sq[bound]) {
                    sl_state.bath_delta[bound] = b_tilda;
                } else if (std::signbit(b_tilda) == std::signbit(b_delta)) {
                    sl_state.bath_delta[bound] = std::copysign(1.0, b_tilda) *
                                                 std::min(std::abs(b_tilda), std::abs(PostProcessing::nu * b_delta));
                } else {
                    sl_state.bath_delta[bound] = 0.0;
                }
            }

            double delta_sum = 0.0;
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                delta_sum += sl_state.bath_delta[bound];
            }
            if (delta_sum != 0.0) {
                double positive = 0.0;
                double negative = 0.0;
                for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                    positive += std::max(0.0, sl_state.bath_delta[bound]);
                    negative += std::max(0.0, -sl_state.bath_delta[bound]);
                }
                const double theta_positive = std::min(1.0, negative / positive);
                const double theta_negative = std::min(1.0, positive / negative);
                for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                    sl_state.bath_delta[bound] = theta_positive * std::max(0.0, sl_state.bath_delta[bound]) -
                                                 theta_negative * std::max(0.0, -sl_state.bath_delta[bound]);
                }
            }

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                sl_state.bath_at_vrtx[vrtx] = sl_state.bath_at_baryctr;
            }
            sl_state.bath_at_vrtx += sl_state.bath_delta * sl_state.T;
            const double del_bath_norm = norm(sl_state.bath_at_vrtx - sl_state.bath_lin);
            const double bath_norm     = norm(sl_state.bath_lin);
            if (del_bath_norm / bath_norm > 1.0e-12) {
                row(elt.data.state[0].aux, SWE::Auxiliaries::bath) = elt.ProjectLinearToBasis(sl_state.bath_at_vrtx);
                row(sl_state.q_lin, SWE::Variables::ze) += (sl_state.bath_lin - sl_state.bath_at_vrtx);
                row(elt.data.state[0].q, SWE::Variables::ze) =
                    elt.ProjectLinearToBasis(row(sl_state.q_lin, SWE::Variables::ze));
            }
        }
    });
}
}

#endif