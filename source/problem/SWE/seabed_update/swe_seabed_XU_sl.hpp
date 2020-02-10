#ifndef SWE_SEABED_XU_SL_HPP
#define SWE_SEABED_XU_SL_HPP

namespace SWE {
template <typename DiscretizationType>
void XU_seabed_slope_limiter(DiscretizationType& discretization) {
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
                const uint element_1 = bound;
                const uint element_2 = (bound + 1) % elt.data.get_nbound();

                StatVector<double, SWE::n_dimensions> del_bath;
                del_bath[0] = sl_state.bath_at_baryctr_neigh[element_1] - sl_state.bath_at_baryctr;
                del_bath[1] = sl_state.bath_at_baryctr_neigh[element_2] - sl_state.bath_at_baryctr;

                column(sl_state.dbath_r, bound) = sl_state.A_inv[bound] * del_bath;
                sl_state.alpha_r[bound] =
                    sl_state.d_r[bound] / (1 + std::pow(sq_norm(column(sl_state.dbath_r, bound)), 2.0));
            }

            const double alpha_r_sum = std::accumulate(sl_state.alpha_r.begin(), sl_state.alpha_r.end(), 0.0);
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                sl_state.w_r[bound] = sl_state.alpha_r[bound] / alpha_r_sum;
            }

            StatVector<double, SWE::n_dimensions> dbath = sl_state.dbath_r * sl_state.w_r;
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                sl_state.bath_delta[bound] = transpose(dbath) * (sl_state.midpts_coord[bound] - sl_state.baryctr_coord);
            }

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                sl_state.bath_at_vrtx[vrtx] = sl_state.bath_at_baryctr;
            }
            sl_state.bath_at_vrtx += sl_state.bath_delta * sl_state.T;
            row(elt.data.state[0].aux, SWE::Auxiliaries::bath) = elt.ProjectLinearToBasis(sl_state.bath_at_vrtx);
            row(sl_state.q_lin, SWE::Variables::ze) += (sl_state.bath_lin - sl_state.bath_at_vrtx);
            row(elt.data.state[0].q, SWE::Variables::ze) =
                elt.ProjectLinearToBasis(row(sl_state.q_lin, SWE::Variables::ze));
        }
    });
}
}

#endif