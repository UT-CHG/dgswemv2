#ifndef RKDG_SWE_POST_WET_DRY_HPP
#define RKDG_SWE_POST_WET_DRY_HPP

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::wetting_drying_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage + 1];
    auto& wd_state = elt.data.wet_dry_state;

    wd_state.q_lin     = elt.ProjectBasisToLinear(state.q);
    wd_state.q_at_vrtx = elt.ComputeLinearUvrtx(wd_state.q_lin);

    for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
        wd_state.h_at_vrtx[vrtx] = wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) + wd_state.bath_at_vrtx[vrtx];
    }

    double h_avg = std::accumulate(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end(), 0.0) / elt.data.get_nvrtx();

    bool set_wet_element = false;
    bool set_dry_element = false;
    bool check_element   = false;

    if (h_avg <= PostProcessing::h_o + PostProcessing::h_o_threshold) {
        // Say that "low" h_avg is 10% of h_o
        // Check if element has dangerously low h_avg or even negative h_avg
        // In this case we need to bump up the water to at least "low" h_avg
        if (h_avg < 0.1 * PostProcessing::h_o) {
            wd_state.went_completely_dry = true;

            h_avg = 0.1 * PostProcessing::h_o;
        }

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) = h_avg - wd_state.bath_at_vrtx[vrtx];
        }

        set_dry_element = true;
    } else {
        uint n_dry_vrtx = 0;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            if (wd_state.h_at_vrtx[vrtx] <= PostProcessing::h_o + PostProcessing::h_o_threshold) {
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

            wd_state.h_at_vrtx_temp[h_min_vrtx] = PostProcessing::h_o;

            wd_state.h_at_vrtx_temp[h_mid_vrtx] =
                std::max(PostProcessing::h_o,
                         wd_state.h_at_vrtx[h_mid_vrtx] -
                             (wd_state.h_at_vrtx_temp[h_min_vrtx] - wd_state.h_at_vrtx[h_min_vrtx]) / 2);

            wd_state.h_at_vrtx_temp[h_max_vrtx] =
                wd_state.h_at_vrtx[h_max_vrtx] -
                (wd_state.h_at_vrtx_temp[h_min_vrtx] - wd_state.h_at_vrtx[h_min_vrtx]) -
                (wd_state.h_at_vrtx_temp[h_mid_vrtx] - wd_state.h_at_vrtx[h_mid_vrtx]);

            wd_state.h_at_vrtx = wd_state.h_at_vrtx_temp;

            double del_qx = 0;
            double del_qy = 0;

            n_dry_vrtx = 0;
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                if (wd_state.h_at_vrtx[vrtx] <= PostProcessing::h_o + PostProcessing::h_o_threshold) {
                    n_dry_vrtx++;

                    del_qx += wd_state.q_at_vrtx(SWE::Variables::qx, vrtx);
                    del_qy += wd_state.q_at_vrtx(SWE::Variables::qy, vrtx);

                    wd_state.q_at_vrtx(SWE::Variables::qx, vrtx) = 0.0;
                    wd_state.q_at_vrtx(SWE::Variables::qy, vrtx) = 0.0;
                }
            }

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) = wd_state.h_at_vrtx[vrtx] - wd_state.bath_at_vrtx[vrtx];

                if (wd_state.h_at_vrtx[vrtx] > PostProcessing::h_o + PostProcessing::h_o_threshold) {
                    wd_state.q_at_vrtx(SWE::Variables::qx, vrtx) += del_qx / (3 - n_dry_vrtx);
                    wd_state.q_at_vrtx(SWE::Variables::qy, vrtx) += del_qy / (3 - n_dry_vrtx);
                }
            }

            state.q = elt.ProjectLinearToBasis(elt.data.get_ndof(), wd_state.q_at_vrtx);

            if (wd_state.wet) {
                set_wet_element = true;
            } else {
                check_element = true;
            }
        }

        if (check_element) {
            auto h_max_iter = std::max_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
            uint h_max_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_max_iter);

            double ze_h_max_vrtx = wd_state.q_at_vrtx(SWE::Variables::ze, h_max_vrtx);

            if (ze_h_max_vrtx > PostProcessing::h_o - wd_state.bath_min) {
                set_wet_element = true;
            } else {
                set_dry_element = true;
            }
        };
    }

    if (set_dry_element) {
        wd_state.wet = false;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.q_at_vrtx(SWE::Variables::qx, vrtx) = 0.0;
            wd_state.q_at_vrtx(SWE::Variables::qy, vrtx) = 0.0;
        }

        state.q = elt.ProjectLinearToBasis(elt.data.get_ndof(), wd_state.q_at_vrtx);

        set_constant(state.rhs, 0.0);
    } else if (set_wet_element) {
        wd_state.wet = true;
    }
}
}
}

#endif