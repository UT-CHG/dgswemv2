#ifndef SWE_PROC_SLOPE_LIMIT_HPP
#define SWE_PROC_SLOPE_LIMIT_HPP

namespace SWE {
void Problem::slope_limiting_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& sl_state = elt.data.slope_limit_state;

        elt.ComputeUbaryctr(state.ze, sl_state.ze_at_baryctr);
        elt.ComputeUbaryctr(state.qx, sl_state.qx_at_baryctr);
        elt.ComputeUbaryctr(state.qy, sl_state.qy_at_baryctr);

        elt.ComputeUmidpts(state.ze, sl_state.ze_at_midpts);
        elt.ComputeUmidpts(state.qx, sl_state.qx_at_midpts);
        elt.ComputeUmidpts(state.qy, sl_state.qy_at_midpts);
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.ze_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.ze_at_baryctr;
        sl_state_in.qx_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.qx_at_baryctr;
        sl_state_in.qy_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.qy_at_baryctr;

        sl_state_ex.ze_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.ze_at_baryctr;
        sl_state_ex.qx_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.qx_at_baryctr;
        sl_state_ex.qy_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.qy_at_baryctr;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.ze_at_baryctr_neigh[bound.bound_id] = sl_state.ze_at_baryctr;
        sl_state.qx_at_baryctr_neigh[bound.bound_id] = sl_state.qx_at_baryctr;
        sl_state.qy_at_baryctr_neigh[bound.bound_id] = sl_state.qy_at_baryctr;
    });

    mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& sl_state = elt.data.slope_limit_state;

        double u = sl_state.qx_at_baryctr / (sl_state.ze_at_baryctr + sl_state.bath_at_baryctr);
        double v = sl_state.qy_at_baryctr / (sl_state.ze_at_baryctr + sl_state.bath_at_baryctr);
        double c = std::sqrt(Global::g * (sl_state.ze_at_baryctr + sl_state.bath_at_baryctr));

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % 3;

            sl_state.R[0][0] = 1.0;
            sl_state.R[1][0] = u + c * sl_state.surface_normal[bound][GlobalCoord::x];
            sl_state.R[2][0] = v + c * sl_state.surface_normal[bound][GlobalCoord::y];

            sl_state.R[0][1] = 0.0;
            sl_state.R[1][1] = -sl_state.surface_normal[bound][GlobalCoord::y];
            sl_state.R[2][1] = sl_state.surface_normal[bound][GlobalCoord::x];

            sl_state.R[0][2] = 1.0;
            sl_state.R[1][2] = u - c * sl_state.surface_normal[bound][GlobalCoord::x];
            sl_state.R[2][2] = v - c * sl_state.surface_normal[bound][GlobalCoord::y];

            double det = sl_state.R[0][0] * sl_state.R[1][1] * sl_state.R[2][2] +
                         sl_state.R[0][1] * sl_state.R[1][2] * sl_state.R[2][0] +
                         sl_state.R[0][2] * sl_state.R[1][0] * sl_state.R[2][1] -
                         sl_state.R[0][0] * sl_state.R[1][2] * sl_state.R[2][1] -
                         sl_state.R[0][1] * sl_state.R[1][0] * sl_state.R[2][2] -
                         sl_state.R[0][2] * sl_state.R[1][1] * sl_state.R[2][0];

            sl_state.L[0][0] = (sl_state.R[1][1] * sl_state.R[2][2] - sl_state.R[1][2] * sl_state.R[2][1]) / det;
            sl_state.L[1][0] = (sl_state.R[1][2] * sl_state.R[2][0] - sl_state.R[1][0] * sl_state.R[2][2]) / det;
            sl_state.L[2][0] = (sl_state.R[1][0] * sl_state.R[2][1] - sl_state.R[1][1] * sl_state.R[2][0]) / det;

            sl_state.L[0][1] = (sl_state.R[0][2] * sl_state.R[2][1] - sl_state.R[0][1] * sl_state.R[2][2]) / det;
            sl_state.L[1][1] = (sl_state.R[0][0] * sl_state.R[2][2] - sl_state.R[0][2] * sl_state.R[2][0]) / det;
            sl_state.L[2][1] = (sl_state.R[0][1] * sl_state.R[2][0] - sl_state.R[0][0] * sl_state.R[2][1]) / det;

            sl_state.L[0][2] = (sl_state.R[0][1] * sl_state.R[1][2] - sl_state.R[0][2] * sl_state.R[1][1]) / det;
            sl_state.L[1][2] = (sl_state.R[0][2] * sl_state.R[1][0] - sl_state.R[0][0] * sl_state.R[1][2]) / det;
            sl_state.L[2][2] = (sl_state.R[0][0] * sl_state.R[1][1] - sl_state.R[0][1] * sl_state.R[1][0]) / det;

            for (uint var = 0; var < 3; var++) {
                sl_state.w_midpt_char[var] = sl_state.L[var][0] * sl_state.ze_at_midpts[bound] +
                                             sl_state.L[var][1] * sl_state.qx_at_midpts[bound] +
                                             sl_state.L[var][2] * sl_state.qy_at_midpts[bound];

                sl_state.w_baryctr_char[var][0] = sl_state.L[var][0] * sl_state.ze_at_baryctr +
                                                  sl_state.L[var][1] * sl_state.qx_at_baryctr +
                                                  sl_state.L[var][2] * sl_state.qy_at_baryctr;

                sl_state.w_baryctr_char[var][1] = sl_state.L[var][0] * sl_state.ze_at_baryctr_neigh[element_1] +
                                                  sl_state.L[var][1] * sl_state.qx_at_baryctr_neigh[element_1] +
                                                  sl_state.L[var][2] * sl_state.qy_at_baryctr_neigh[element_1];

                sl_state.w_baryctr_char[var][2] = sl_state.L[var][0] * sl_state.ze_at_baryctr_neigh[element_2] +
                                                  sl_state.L[var][1] * sl_state.qx_at_baryctr_neigh[element_2] +
                                                  sl_state.L[var][2] * sl_state.qy_at_baryctr_neigh[element_2];
            }

            double w_tilda;
            double w_delta;

            double M = 50;
            double nu = 1.5;

            for (uint var = 0; var < 3; var++) {
                w_tilda = sl_state.w_midpt_char[var] - sl_state.w_baryctr_char[var][0];

                w_delta =
                    sl_state.alpha_1[bound] * (sl_state.w_baryctr_char[var][1] - sl_state.w_baryctr_char[var][0]) +
                    sl_state.alpha_2[bound] * (sl_state.w_baryctr_char[var][2] - sl_state.w_baryctr_char[var][0]);

                // TVB modified minmod
                if (std::abs(w_tilda) <= M * sl_state.r_sq[bound]) {
                    sl_state.delta_char[var] = w_tilda;
                } else if (std::signbit(w_tilda) == std::signbit(w_delta)) {
                    sl_state.delta_char[var] =
                        std::copysign(1.0, w_tilda) * std::min(std::abs(w_tilda), std::abs(nu * w_delta));
                } else {
                    sl_state.delta_char[var] = 0.0;
                }
            }

            for (uint var = 0; var < 3; var++) {
                sl_state.delta[var][bound] = sl_state.R[var][0] * sl_state.delta_char[0] +
                                             sl_state.R[var][1] * sl_state.delta_char[1] +
                                             sl_state.R[var][2] * sl_state.delta_char[2];
            }
        }

        for (uint var = 0; var < 3; var++) {
            double delta_sum = 0.0;

            for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
                delta_sum += sl_state.delta[var][bound];
            }

            if (delta_sum != 0.0) {
                double positive = 0.0;
                double negative = 0.0;

                for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
                    positive += std::max(0.0, sl_state.delta[var][bound]);
                    negative += std::max(0.0, -sl_state.delta[var][bound]);
                }

                double theta_positive = std::min(1.0, negative / positive);
                double theta_negative = std::min(1.0, positive / negative);

                for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
                    sl_state.delta[var][bound] = theta_positive * std::max(0.0, sl_state.delta[var][bound]) -
                                                 theta_negative * std::max(0.0, -sl_state.delta[var][bound]);
                }
            }
        }

        Array2D<double> T{{-1.0, 1.0, 1.0}, {1.0, -1.0, 1.0}, {1.0, 1.0, -1.0}};

        for (uint vrtx = 0; vrtx < 3; vrtx++) {
            sl_state.ze_at_vrtx[vrtx] = sl_state.ze_at_baryctr;
            sl_state.qx_at_vrtx[vrtx] = sl_state.qx_at_baryctr;
            sl_state.qy_at_vrtx[vrtx] = sl_state.qy_at_baryctr;

            for (uint bound = 0; elt.data.get_nbound(); bound++) {
                sl_state.ze_at_vrtx[vrtx] += T[vrtx][bound] * sl_state.delta[0][bound];
                sl_state.qx_at_vrtx[vrtx] += T[vrtx][bound] * sl_state.delta[1][bound];
                sl_state.qy_at_vrtx[vrtx] += T[vrtx][bound] * sl_state.delta[2][bound];
            }
        }

        state.ze = elt.L2Projection(sl_state.ze_at_vrtx);
        state.qx = elt.L2Projection(sl_state.qx_at_vrtx);
        state.qy = elt.L2Projection(sl_state.qy_at_vrtx);
    });
}
}

#endif