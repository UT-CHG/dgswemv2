#ifndef SWE_XU_SL_HPP
#define SWE_XU_SL_HPP

// Implementation of XU etal. slope limiter

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

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = bound;
            const uint element_2 = (bound + 1) % elt.data.get_nbound();

            StatMatrix<double, SWE::n_dimensions, SWE::n_variables> del_q;
            row(del_q, 0) = transpose(sl_state.q_at_baryctr_neigh[element_1] - sl_state.q_at_baryctr);
            row(del_q, 1) = transpose(sl_state.q_at_baryctr_neigh[element_2] - sl_state.q_at_baryctr);

            sl_state.dq_r[bound]    = sl_state.A_inv[bound] * del_q;
            sl_state.alpha_r[bound] = sl_state.d_r[bound] / (1 + std::pow(sq_norm(sl_state.dq_r[bound]), 2.0));
        }

        const double alpha_r_sum = std::accumulate(sl_state.alpha_r.begin(), sl_state.alpha_r.end(), 0.0);
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            sl_state.w_r[bound] = sl_state.alpha_r[bound] / alpha_r_sum;
        }

        StatMatrix<double, SWE::n_dimensions, SWE::n_variables> dq;
        set_constant(dq, 0.0);
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            dq += sl_state.dq_r[bound] * sl_state.w_r[bound];
        }
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(sl_state.delta, bound) = transpose(dq) * (sl_state.midpts_coord[bound] - sl_state.baryctr_coord);
        }

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
            column(sl_state.q_at_vrtx, vrtx) = sl_state.q_at_baryctr;
        }
        sl_state.q_at_vrtx += sl_state.delta * sl_state.T;
        state.q = elt.ProjectLinearToBasis(sl_state.q_at_vrtx);
    }
}
}

#endif