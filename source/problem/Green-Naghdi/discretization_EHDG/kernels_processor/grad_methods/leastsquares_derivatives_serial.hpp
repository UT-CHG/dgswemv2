#ifndef LEASTSQUARES_DERIVATIVES_SERIAL_HPP
#define LEASTSQUARES_DERIVATIVES_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dze_ls(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_du_ls(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_ddu_ls(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);

void Problem::compute_derivatives_serial(ProblemDiscretizationType& discretization,
                                         ProblemGlobalDataType& global_data,
                                         const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);
        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        derivative.ze_lin        = elt.ProjectBasisToLinear(row(state.q, SWE::Variables::ze));
        derivative.ze_at_baryctr = elt.ComputeLinearUbaryctr(derivative.ze_lin);
        derivative.ze_at_midpts  = elt.ComputeLinearUmidpts(derivative.ze_lin);
    });
    compute_dze_ls(discretization, stepper);
    reconstruct_dze(discretization, global_data, stepper);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        row(internal.u_at_gp, GlobalCoord::x) =
            vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
        row(internal.u_at_gp, GlobalCoord::y) =
            vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

        derivative.u_lin        = elt.ProjectBasisToLinear(elt.L2Projection(internal.u_at_gp));
        derivative.u_at_baryctr = elt.ComputeLinearUbaryctr(derivative.u_lin);
        derivative.u_at_midpts  = elt.ComputeLinearUmidpts(derivative.u_lin);
    });
    compute_du_ls(discretization, stepper);
    reconstruct_du(discretization, global_data, stepper);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& derivative         = elt.data.derivative;
        derivative.du_at_baryctr = elt.ComputeLinearUbaryctr(derivative.du_lin);
        derivative.du_at_midpts  = elt.ComputeLinearUmidpts(derivative.du_lin);
    });
    compute_ddu_ls(discretization, stepper);
    reconstruct_ddu(discretization, global_data, stepper);
}

template <typename ProblemDiscretizationType>
void compute_dze_ls(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage    = stepper.GetStage();
        auto& state_in      = intface.data_in.state[stage];
        auto& state_ex      = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);
        row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);

        derivative_in.ze_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.ze_at_baryctr;
        derivative_ex.ze_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.ze_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        boundary.q_at_gp = bound.ComputeUgp(state.q);
        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        derivative.ze_at_baryctr_neigh[bound.bound_id] = derivative.ze_at_midpts[bound.bound_id];
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        const uint stage = stepper.GetStage();
        auto& state      = dbound.data.state[stage];
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        boundary.q_at_gp = dbound.ComputeUgp(state.q);
        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        std::vector<double> message(3);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        derivative.ze_at_baryctr_neigh[dbound.bound_id]                = message[0];
        derivative.u_at_baryctr_neigh[dbound.bound_id][GlobalCoord::x] = message[1];
        derivative.u_at_baryctr_neigh[dbound.bound_id][GlobalCoord::y] = message[2];
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            derivative.ze_at_midpts[bound] = derivative.ze_at_baryctr_neigh[bound] - derivative.ze_at_baryctr;
        }
        derivative.dze_at_baryctr = derivative.P * transpose(derivative.ze_at_midpts);
    });
}

template <typename ProblemDiscretizationType>
void compute_du_ls(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& derivative_in                                   = intface.data_in.derivative;
        auto& derivative_ex                                   = intface.data_ex.derivative;
        derivative_in.u_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.u_at_baryctr;
        derivative_ex.u_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.u_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& derivative                              = bound.data.derivative;
        derivative.u_at_baryctr_neigh[bound.bound_id] = column(derivative.u_at_midpts, bound.bound_id);
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.u_at_midpts, bound) = derivative.u_at_baryctr_neigh[bound] - derivative.u_at_baryctr;
        }
        derivative.du_at_baryctr = flatten<double, GN::n_dimensions, GN::n_dimensions, SO::ColumnMajor>(
            derivative.P * transpose(derivative.u_at_midpts));
    });
}

template <typename ProblemDiscretizationType>
void compute_ddu_ls(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& derivative_in                                    = intface.data_in.derivative;
        auto& derivative_ex                                    = intface.data_ex.derivative;
        derivative_in.du_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.du_at_baryctr;
        derivative_ex.du_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.du_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& derivative                               = bound.data.derivative;
        derivative.du_at_baryctr_neigh[bound.bound_id] = column(derivative.du_at_midpts, bound.bound_id);
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
        std::vector<double> message(GN::n_du_terms + ngp);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        for (uint du = 0; du < GN::n_du_terms; ++du) {
            derivative.du_at_baryctr_neigh[dbound.bound_id][du] = message[du];
        }
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            boundary.aux_at_gp(SWE::Auxiliaries::h, gp) =
                (boundary.aux_at_gp(SWE::Auxiliaries::h, gp) + message[GN::n_du_terms + gp_ex]) / 2.0;
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.du_at_midpts, bound) = derivative.du_at_baryctr_neigh[bound] - derivative.du_at_baryctr;
        }
        derivative.ddu_at_baryctr = flatten<double, GN::n_dimensions, GN::n_du_terms, SO::ColumnMajor>(
            derivative.P * transpose(derivative.du_at_midpts));
    });
}
}
}

#endif
