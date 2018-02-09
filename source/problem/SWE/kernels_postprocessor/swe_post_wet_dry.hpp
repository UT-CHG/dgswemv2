#ifndef SWE_POST_WET_DRY_HPP
#define SWE_POST_WET_DRY_HPP

namespace SWE {
template <typename ElementType>
void Problem::wetting_drying_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage + 1];
    auto& wd_state = elt.data.wet_dry_state;
    auto& internal = elt.data.internal;

    elt.ProjectBasisToLinear(state.ze, wd_state.ze_lin);

    elt.ComputeLinearUvrtx(wd_state.ze_lin, wd_state.ze_at_vrtx);

    for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
        wd_state.h_at_vrtx[vrtx] = wd_state.ze_at_vrtx[vrtx] + wd_state.bath_at_vrtx[vrtx];
    }

    double h_avg = std::accumulate(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end(), 0.0) / elt.data.get_nvrtx();

    bool set_wet_element = false;
    bool set_dry_element = false;
    bool check_element = false;

    if (h_avg <= Global::h_o) {
        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.ze_at_vrtx[vrtx] = h_avg - wd_state.bath_at_vrtx[vrtx];
        }
    
        elt.ProjectLinearToBasis(wd_state.ze_at_vrtx, state.ze);

        set_dry_element = true;
    } else {
        uint n_dry_vrtx = 0;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            if (wd_state.h_at_vrtx[vrtx] <= Global::h_o) {
                n_dry_vrtx++;
            }
        }

        if (n_dry_vrtx == 0) {
            if (wd_state.wet) {
                set_wet_element = true;
            } else {
                check_element = true;
            }
        } else {
            auto h_min_iter = std::min_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
            uint h_min_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_min_iter);

            auto h_max_iter = std::max_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
            uint h_max_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_max_iter);

            uint h_mid_vrtx = 3 - h_max_vrtx - h_min_vrtx;

            wd_state.h_at_vrtx_temp[h_min_vrtx] = Global::h_o;

            wd_state.h_at_vrtx_temp[h_mid_vrtx] =
                std::max(Global::h_o,
                         wd_state.h_at_vrtx[h_mid_vrtx] -
                             (wd_state.h_at_vrtx_temp[h_min_vrtx] - wd_state.h_at_vrtx[h_min_vrtx]) / 2);

            wd_state.h_at_vrtx_temp[h_max_vrtx] =
                wd_state.h_at_vrtx[h_max_vrtx] -
                (wd_state.h_at_vrtx_temp[h_min_vrtx] - wd_state.h_at_vrtx[h_min_vrtx]) -
                (wd_state.h_at_vrtx_temp[h_mid_vrtx] - wd_state.h_at_vrtx[h_mid_vrtx]);

            wd_state.h_at_vrtx = wd_state.h_at_vrtx_temp;

            elt.ProjectBasisToLinear(state.qx, wd_state.qx_lin);
            elt.ProjectBasisToLinear(state.qy, wd_state.qy_lin);

            elt.ComputeLinearUvrtx(wd_state.qx_lin, wd_state.qx_at_vrtx);
            elt.ComputeLinearUvrtx(wd_state.qy_lin, wd_state.qy_at_vrtx);

            double del_qx = 0;
            double del_qy = 0;

            n_dry_vrtx = 0;
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                if (wd_state.h_at_vrtx[vrtx] <= Global::h_o) {
                    n_dry_vrtx++;

                    del_qx += wd_state.qx_at_vrtx[vrtx];
                    del_qy += wd_state.qy_at_vrtx[vrtx];

                    wd_state.qx_at_vrtx[vrtx] = 0;
                    wd_state.qy_at_vrtx[vrtx] = 0;
                }
            }

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                wd_state.ze_at_vrtx[vrtx] = wd_state.h_at_vrtx[vrtx] - wd_state.bath_at_vrtx[vrtx];

                if (wd_state.h_at_vrtx[vrtx] > Global::h_o) {
                    wd_state.qx_at_vrtx[vrtx] += del_qx / (3 - n_dry_vrtx);
                    wd_state.qy_at_vrtx[vrtx] += del_qy / (3 - n_dry_vrtx);
                }
            }

            elt.ProjectLinearToBasis(wd_state.ze_at_vrtx, state.ze);
            elt.ProjectLinearToBasis(wd_state.qx_at_vrtx, state.qx);
            elt.ProjectLinearToBasis(wd_state.qy_at_vrtx, state.qy);

            if (wd_state.wet) {
                set_wet_element = true;
            } else {
                check_element = true;
            }
        }

        if (check_element) {
            auto h_max_iter = std::max_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
            uint h_max_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_max_iter);

            double ze_h_max_vrtx = wd_state.ze_at_vrtx[h_max_vrtx];

            if (ze_h_max_vrtx > Global::h_o - wd_state.bath_min) {
                set_wet_element = true;
            } else {
                set_dry_element = true;
            }
        };
    }

    if (set_dry_element) {
        wd_state.wet = false;

        std::fill(state.qx.begin(), state.qx.end(), 0.0);
        std::fill(state.qy.begin(), state.qy.end(), 0.0);

        std::fill(state.rhs_ze.begin(), state.rhs_ze.end(), 0.0);
        std::fill(state.rhs_qx.begin(), state.rhs_qx.end(), 0.0);
        std::fill(state.rhs_qy.begin(), state.rhs_qy.end(), 0.0);
    } else if (set_wet_element) {
        wd_state.wet = true;
    }

    elt.ComputeUgp(state.ze, internal.ze_at_gp);

    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.h_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];
    }

    wd_state.water_volume = elt.Integration(internal.h_at_gp);
}
}

#endif