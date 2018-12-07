#ifndef SWE_CS_SL_HPP
#define SWE_CS_SL_HPP

// Implementation of Cockburn-Shu slope limiter
#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"

namespace SWE {
template <typename StepperType, typename ElementType>
void slope_limiting_prepare_element_kernel(const StepperType& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;
    auto& sl_state = elt.data.slope_limit_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        sl_state.q_lin        = elt.ProjectBasisToLinear(state.q);
        sl_state.q_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.q_lin);
        sl_state.q_at_midpts  = elt.ComputeLinearUmidpts(sl_state.q_lin);
    }
}

template <typename StepperType, typename InterfaceType>
void slope_limiting_prepare_interface_kernel(const StepperType& stepper, InterfaceType& intface) {
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

template <typename StepperType, typename BoundaryType>
void slope_limiting_prepare_boundary_kernel(const StepperType& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;
    auto& sl_state = bound.data.slope_limit_state;

    sl_state.wet_neigh[bound.bound_id] = wd_state.wet;

    if (wd_state.wet) {
        sl_state.q_at_baryctr_neigh[bound.bound_id] = sl_state.q_at_baryctr;
    }
}

template <typename StepperType, typename DistributedBoundaryType>
void slope_limiting_distributed_boundary_send_kernel(const StepperType& stepper,
                                                     DistributedBoundaryType& dbound,
                                                     uint comm_type) {
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
    dbound.boundary_condition.exchanger.SetToSendBuffer(comm_type, message);
}

template <typename StepperType, typename DistributedBoundaryType>
void slope_limiting_prepare_distributed_boundary_kernel(const StepperType& stepper,
                                                        DistributedBoundaryType& dbound,
                                                        uint comm_type) {
    auto& sl_state = dbound.data.slope_limit_state;

    std::vector<double> message;

    message.resize(1 + SWE::n_variables);

    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(comm_type, message);

    sl_state.wet_neigh[dbound.bound_id] = message[0];

    for (uint var = 0; var < SWE::n_variables; ++var) {
        sl_state.q_at_baryctr_neigh[dbound.bound_id][var] = message[1 + var];
    }
}

template <typename StepperType, typename ElementType>
void slope_limiting_kernel(const StepperType& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;
    auto& sl_state = elt.data.slope_limit_state;

    if (wd_state.wet &&
        std::find(sl_state.wet_neigh.begin(), sl_state.wet_neigh.end(), false) == sl_state.wet_neigh.end()) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        StatMatrix<double, SWE::n_variables, SWE::n_variables> R;
        StatMatrix<double, SWE::n_variables, SWE::n_variables> invR;

        StatVector<double, SWE::n_variables> w_midpt_char;
        StatVector<double, SWE::n_variables> w_baryctr_char_0;
        StatVector<double, SWE::n_variables> w_baryctr_char_1;
        StatVector<double, SWE::n_variables> w_baryctr_char_2;

        StatVector<double, SWE::n_variables> w_tilda;
        StatVector<double, SWE::n_variables> w_delta;
        StatVector<double, SWE::n_variables> delta_char;

        double h = sl_state.q_at_baryctr[SWE::Variables::ze] + sl_state.bath_at_baryctr;
        double u = sl_state.q_at_baryctr[SWE::Variables::qx] / h;
        double v = sl_state.q_at_baryctr[SWE::Variables::qy] / h;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % 3;

            double nx = sl_state.surface_normal[bound][GlobalCoord::x];
            double ny = sl_state.surface_normal[bound][GlobalCoord::y];

            R    = SWE::R(h, u, v, nx, ny);
            invR = SWE::invR(h, u, v, nx, ny);

            w_midpt_char     = invR * column(sl_state.q_at_midpts, bound);
            w_baryctr_char_0 = invR * sl_state.q_at_baryctr;
            w_baryctr_char_1 = invR * sl_state.q_at_baryctr_neigh[element_1];
            w_baryctr_char_2 = invR * sl_state.q_at_baryctr_neigh[element_2];

            w_tilda = w_midpt_char - w_baryctr_char_0;

            w_delta = sl_state.alpha_1[bound] * (w_baryctr_char_1 - w_baryctr_char_0) +
                      sl_state.alpha_2[bound] * (w_baryctr_char_2 - w_baryctr_char_0);

            for (uint var = 0; var < SWE::n_variables; ++var) {
                // TVB modified minmod
                if (std::abs(w_tilda[var]) <= PostProcessing::M * sl_state.r_sq[bound]) {
                    delta_char[var] = w_tilda[var];
                } else if (std::signbit(w_tilda[var]) == std::signbit(w_delta[var])) {
                    delta_char[var] = std::copysign(1.0, w_tilda[var]) *
                                      std::min(std::abs(w_tilda[var]), std::abs(PostProcessing::nu * w_delta[var]));
                } else {
                    delta_char[var] = 0.0;
                }
            }

            column(sl_state.delta, bound) = R * delta_char;
        }

        for (uint var = 0; var < SWE::n_variables; ++var) {
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

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
            column(sl_state.q_at_vrtx, vrtx) = sl_state.q_at_baryctr;
        }

        StatMatrix<double, SWE::n_variables, SWE::n_variables> T;
        set_constant(T, 1.0);
        T(0, 0) = -1.0;
        T(1, 1) = -1.0;
        T(2, 2) = -1.0;

        sl_state.q_at_vrtx += sl_state.delta * T;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            double del_q_norm = norm(row(sl_state.q_at_vrtx, var) - row(sl_state.q_lin, var));
            double q_norm     = norm(row(sl_state.q_lin, var));

            if (del_q_norm / q_norm > 1.0e-12) {
                row(state.q, var) = elt.ProjectLinearToBasis(elt.data.get_ndof(), row(sl_state.q_at_vrtx, var));
            }
        }
    }
}
}

#endif