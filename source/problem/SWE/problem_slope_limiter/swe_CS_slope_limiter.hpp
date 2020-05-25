#ifndef SWE_CS_SL_HPP
#define SWE_CS_SL_HPP

// Implementation of Cockburn-Shu slope limiter
#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"

namespace SWE {
template <typename StepperType, typename ElementType>
void slope_limiting_prepare_element_kernel(const StepperType& stepper, ElementType& elt) {
    auto& state    = elt.data.state[stepper.GetStage()];
    auto& wd_state = elt.data.wet_dry_state;
    auto& sl_state = elt.data.slope_limit_state;
    if (wd_state.wet) {
        sl_state.q_lin        = elt.ProjectBasisToLinear(state.q);
        sl_state.q_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.q_lin);
        sl_state.q_at_vrtx    = sl_state.q_lin;
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
        sl_state.q_at_baryctr_neigh[bound.bound_id] = column(sl_state.q_at_midpts, bound.bound_id);
    }
}

template <typename StepperType, typename DistributedBoundaryType>
void slope_limiting_distributed_boundary_send_kernel(const StepperType& stepper,
                                                     DistributedBoundaryType& dbound,
                                                     uint comm_type) {
    auto& state    = dbound.data.state[stepper.GetStage()];
    auto& boundary = dbound.data.boundary[dbound.bound_id];
    auto& wd_state = dbound.data.wet_dry_state;
    auto& sl_state = dbound.data.slope_limit_state;

    const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    std::vector<double> message(1 + SWE::n_variables + ngp);
    message[0] = (double)wd_state.wet;
    if (wd_state.wet) {
        boundary.q_at_gp = dbound.ComputeUgp(state.q);
        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);
        for (uint var = 0; var < SWE::n_variables; ++var) {
            message[1 + var] = sl_state.q_at_baryctr[var];
        }
        for (uint gp = 0; gp < ngp; ++gp) {
            message[1 + SWE::n_variables + gp] = boundary.aux_at_gp(SWE::Auxiliaries::h, gp);
        }
    }
    dbound.boundary_condition.exchanger.SetToSendBuffer(comm_type, message);
}

template <typename StepperType, typename DistributedBoundaryType>
void slope_limiting_prepare_distributed_boundary_kernel(const StepperType& stepper,
                                                        DistributedBoundaryType& dbound,
                                                        uint comm_type) {
    auto& sl_state = dbound.data.slope_limit_state;

    std::vector<double> message(1 + SWE::n_variables);
    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(comm_type, message);
    sl_state.wet_neigh[dbound.bound_id] = (bool)message[0];
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
        auto& state = elt.data.state[stepper.GetStage()];

        StatMatrix<double, SWE::n_variables, SWE::n_variables> R;
        StatMatrix<double, SWE::n_variables, SWE::n_variables> invR;

        StatVector<double, SWE::n_variables> w_midpt_char;
        StatVector<double, SWE::n_variables> w_baryctr_char_0;
        StatVector<double, SWE::n_variables> w_baryctr_char_1;
        StatVector<double, SWE::n_variables> w_baryctr_char_2;

        StatVector<double, SWE::n_variables> w_tilda;
        StatVector<double, SWE::n_variables> w_delta;
        StatVector<double, SWE::n_variables> delta_char;

        const double h = sl_state.q_at_baryctr[SWE::Variables::ze] + sl_state.bath_at_baryctr;
        const double u = sl_state.q_at_baryctr[SWE::Variables::qx] / h;
        const double v = sl_state.q_at_baryctr[SWE::Variables::qy] / h;
        const double c = sl_state.q_at_baryctr[SWE::Variables::hc] / h;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = sl_state.a_elem[2 * bound];
            const uint element_2 = sl_state.a_elem[2 * bound + 1];
            const double mx      = sl_state.median[bound][GlobalCoord::x];
            const double my      = sl_state.median[bound][GlobalCoord::y];

            const SWE::parameters param{h, u, v, c, mx, my};
            R    = SWE::R(param);
            invR = SWE::invR(param);

            w_midpt_char     = invR * column(sl_state.q_at_midpts, bound);
            w_baryctr_char_0 = invR * sl_state.q_at_baryctr;
            w_baryctr_char_1 = invR * sl_state.q_at_baryctr_neigh[element_1];
            w_baryctr_char_2 = invR * sl_state.q_at_baryctr_neigh[element_2];

            w_tilda = w_midpt_char - w_baryctr_char_0;
            w_delta = sl_state.alpha[bound][0] * (w_baryctr_char_1 - w_baryctr_char_0) +
                      sl_state.alpha[bound][1] * (w_baryctr_char_2 - w_baryctr_char_0);
            for (uint var = 0; var < SWE::n_variables; ++var) {
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
        sl_state.q_at_vrtx += sl_state.delta * sl_state.T;
        for (uint var = 0; var < SWE::n_variables; ++var) {
            double del_q_norm = norm(row(sl_state.q_at_vrtx, var) - row(sl_state.q_lin, var));
            double q_norm     = norm(row(sl_state.q_lin, var));
            if (del_q_norm / q_norm > 1.0e-12) {
                row(state.q, var) = elt.ProjectLinearToBasis(row(sl_state.q_at_vrtx, var));
            }
        }
    }
}
}

#endif
