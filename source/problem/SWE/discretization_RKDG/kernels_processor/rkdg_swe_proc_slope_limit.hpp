#ifndef RKDG_SWE_PROC_SLOPE_LIMIT_HPP
#define RKDG_SWE_PROC_SLOPE_LIMIT_HPP

#include "utilities/almost_equal.hpp"

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::slope_limiting_prepare_element_kernel(const RKStepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;
    auto& sl_state = elt.data.slope_limit_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        sl_state.q_lin        = elt.ProjectBasisToLinear(state.q);
        sl_state.q_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.q_lin);
        sl_state.q_at_midpts  = elt.ComputeLinearUmidpts(sl_state.q_lin);
    }
}

template <typename InterfaceType>
void Problem::slope_limiting_prepare_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    auto& wd_state_in = intface.data_in.wet_dry_state;
    auto& wd_state_ex = intface.data_ex.wet_dry_state;

    auto& sl_state_in = intface.data_in.slope_limit_state;
    auto& sl_state_ex = intface.data_ex.slope_limit_state;

    sl_state_in.wet_neigh[intface.bound_id_in] = wd_state_ex.wet;
    sl_state_ex.wet_neigh[intface.bound_id_ex] = wd_state_in.wet;

    if (wd_state_in.wet && wd_state_ex.wet) {
        sl_state_in.q_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.q_at_baryctr;
        sl_state_ex.q_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.q_at_baryctr;
    }
}

template <typename BoundaryType>
void Problem::slope_limiting_prepare_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;
    auto& sl_state = bound.data.slope_limit_state;

    sl_state.wet_neigh[bound.bound_id] = wd_state.wet;

    if (wd_state.wet) {
        sl_state.q_at_baryctr_neigh[bound.bound_id] = sl_state.q_at_baryctr;
    }
}

template <typename DistributedBoundaryType>
void Problem::slope_limiting_distributed_boundary_send_kernel(const RKStepper& stepper,
                                                              DistributedBoundaryType& dbound) {
    auto& wd_state = dbound.data.wet_dry_state;

    // Construct message to exterior state
    std::vector<double> message;

    message.reserve(1 + SWE::n_variables);

    message.push_back(wd_state.wet);

    if (wd_state.wet) {
        auto& sl_state = dbound.data.slope_limit_state;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            message.push_back(sl_state.q_at_baryctr[var]);
        }
    }

    // Set message to send buffer
    dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::baryctr_state, message);
}

template <typename DistributedBoundaryType>
void Problem::slope_limiting_prepare_distributed_boundary_kernel(const RKStepper& stepper,
                                                                 DistributedBoundaryType& dbound) {
    auto& sl_state = dbound.data.slope_limit_state;

    std::vector<double> message;

    message.resize(1 + SWE::n_variables);

    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::baryctr_state, message);

    sl_state.wet_neigh[dbound.bound_id] = message[0];

    for (uint var = 0; var < SWE::n_variables; ++var) {
        sl_state.q_at_baryctr_neigh[dbound.bound_id][var] = message[1 + var];
    }
}

template <typename ElementType>
void Problem::slope_limiting_kernel(const RKStepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;
    auto& sl_state = elt.data.slope_limit_state;

    if (wd_state.wet &&
        std::find(sl_state.wet_neigh.begin(), sl_state.wet_neigh.end(), false) == sl_state.wet_neigh.end()) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        double u = sl_state.q_at_baryctr[SWE::Variables::qx] /
                   (sl_state.q_at_baryctr[SWE::Variables::ze] + sl_state.bath_at_baryctr);
        double v = sl_state.q_at_baryctr[SWE::Variables::qy] /
                   (sl_state.q_at_baryctr[SWE::Variables::ze] + sl_state.bath_at_baryctr);
        double c = std::sqrt(Global::g * (sl_state.q_at_baryctr[SWE::Variables::ze] + sl_state.bath_at_baryctr));

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % 3;

            sl_state.R(0, 0) = 1.0;
            sl_state.R(1, 0) = u + c * sl_state.surface_normal[bound][GlobalCoord::x];
            sl_state.R(2, 0) = v + c * sl_state.surface_normal[bound][GlobalCoord::y];

            sl_state.R(0, 1) = 0.0;
            sl_state.R(1, 1) = -sl_state.surface_normal[bound][GlobalCoord::y];
            sl_state.R(2, 1) = sl_state.surface_normal[bound][GlobalCoord::x];

            sl_state.R(0, 2) = 1.0;
            sl_state.R(1, 2) = u - c * sl_state.surface_normal[bound][GlobalCoord::x];
            sl_state.R(2, 2) = v - c * sl_state.surface_normal[bound][GlobalCoord::y];

            sl_state.L = inverse(sl_state.R);

            sl_state.w_midpt_char = sl_state.L * column(sl_state.q_at_midpts, bound);

            column(sl_state.w_baryctr_char, 0) = sl_state.L * sl_state.q_at_baryctr;
            column(sl_state.w_baryctr_char, 1) = sl_state.L * sl_state.q_at_baryctr_neigh[element_1];
            column(sl_state.w_baryctr_char, 2) = sl_state.L * sl_state.q_at_baryctr_neigh[element_2];

            double w_tilda;
            double w_delta;

            for (uint var = 0; var < 3; var++) {
                w_tilda = sl_state.w_midpt_char[var] - sl_state.w_baryctr_char(var, 0);

                w_delta =
                    sl_state.alpha_1[bound] * (sl_state.w_baryctr_char(var, 1) - sl_state.w_baryctr_char(var, 0)) +
                    sl_state.alpha_2[bound] * (sl_state.w_baryctr_char(var, 2) - sl_state.w_baryctr_char(var, 0));

                // TVB modified minmod
                if (std::abs(w_tilda) <= PostProcessing::M * sl_state.r_sq[bound]) {
                    sl_state.delta_char[var] = w_tilda;
                } else if (std::signbit(w_tilda) == std::signbit(w_delta)) {
                    sl_state.delta_char[var] = std::copysign(1.0, w_tilda) *
                                               std::min(std::abs(w_tilda), std::abs(PostProcessing::nu * w_delta));
                } else {
                    sl_state.delta_char[var] = 0.0;
                }
            }

            column(sl_state.delta, bound) = sl_state.R * sl_state.delta_char;
        }

        for (uint var = 0; var < 3; var++) {
            double delta_sum = 0.0;

            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                delta_sum += sl_state.delta(var, bound);
            }

            if (delta_sum != 0.0) {
                double positive = 0.0;
                double negative = 0.0;

                for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                    positive += std::max(0.0, sl_state.delta(var, bound));
                    negative += std::max(0.0, -sl_state.delta(var, bound));
                }

                double theta_positive = std::min(1.0, negative / positive);
                double theta_negative = std::min(1.0, positive / negative);

                for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                    sl_state.delta(var, bound) = theta_positive * std::max(0.0, sl_state.delta(var, bound)) -
                                                 theta_negative * std::max(0.0, -sl_state.delta(var, bound));
                }
            }
        }

        StatMatrix<double, SWE::n_variables, SWE::n_variables> T;

        set_constant(T, 1.0);

        T(0, 0) = -1.0;
        T(1, 1) = -1.0;
        T(2, 2) = -1.0;

        for (uint vrtx = 0; vrtx < 3; ++vrtx) {
            column(sl_state.q_at_vrtx, vrtx) = sl_state.q_at_baryctr;

            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                sl_state.q_at_vrtx(SWE::Variables::ze, vrtx) += T(vrtx, bound) * sl_state.delta(0, bound);
                sl_state.q_at_vrtx(SWE::Variables::qx, vrtx) += T(vrtx, bound) * sl_state.delta(1, bound);
                sl_state.q_at_vrtx(SWE::Variables::qy, vrtx) += T(vrtx, bound) * sl_state.delta(2, bound);
            }
        }

        bool limit = false;

        for (uint vrtx = 0; vrtx < 3; ++vrtx) {
            double del_q_norm = norm(column(sl_state.q_at_vrtx, vrtx) - column(sl_state.q_lin, vrtx));
            double q_norm     = norm(column(sl_state.q_lin, vrtx));

            if (del_q_norm / q_norm > 1.0e-6) {
                limit = true;
            }
        }

        if (limit) {
            state.q = elt.ProjectLinearToBasis(elt.data.get_ndof(), sl_state.q_at_vrtx);
        }
    }
}
}
}

#endif