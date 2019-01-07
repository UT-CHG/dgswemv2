#ifndef IHDG_SWE_BC_FUNCTION_HPP
#define IHDG_SWE_BC_FUNCTION_HPP

namespace SWE {
namespace IHDG {
namespace BC {
class Function {
  private:
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> Aplus;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dqy;

    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> Aminus;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dqy;

  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename EdgeBoundaryType>
    void ComputeInitTrace(EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename BoundaryType>
void Function::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);

    this->Aplus.resize(ngp);
    this->dAplus_dze.resize(ngp);
    this->dAplus_dqx.resize(ngp);
    this->dAplus_dqy.resize(ngp);

    this->Aminus.resize(ngp);
    this->dAminus_dze.resize(ngp);
    this->dAminus_dqx.resize(ngp);
    this->dAminus_dqy.resize(ngp);
}

template <typename EdgeBoundaryType>
void Function::ComputeInitTrace(EdgeBoundaryType& edge_bound) {
    auto& bound = edge_bound.boundary;

    auto& state    = bound.data.state[0];
    auto& boundary = bound.data.boundary[bound.bound_id];

    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    boundary.q_at_gp = bound.ComputeUgp(state.q);

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        get_Aplus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aplus);
        get_dAplus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dze);
        get_dAplus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqx);
        get_dAplus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqy);

        get_Aminus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aminus);
        get_dAminus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dze);
        get_dAminus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqx);
        get_dAminus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqy);

        HybMatrix<double, SWE::n_variables> q_ex = edge_bound.boundary.ComputeFgp([](Point<2>& pt) {
            double t = 0.0;

            double ze = 0.0;
            double qx = 0.0;
            double qy = 0.0;

            Utilities::ignore(ze, qx, qy);

            if (t <= 3.0) {
                ze = cos(PI * t) - 1.0;
            } else {
                ze = -2.0;
            }

            // StatVector<double, SWE::n_variables> q{ze, qx, qy};
            StatVector<double, SWE::n_variables> q(SWE::ic_q(t, pt));

            return q;
        });

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            auto q     = column(boundary.q_at_gp, gp);
            auto q_hat = column(edge_internal.q_hat_at_gp, gp);
            auto q_inf = column(q_ex, gp);

            StatMatrix<double, SWE::n_variables, SWE::n_variables> dB_dq_hat = this->Aminus[gp] - this->Aplus[gp];

            column(dB_dq_hat, SWE::Variables::ze) +=
                this->dAplus_dze[gp] * (q - q_hat) - this->dAminus_dze[gp] * (q_inf - q_hat);
            column(dB_dq_hat, SWE::Variables::qx) +=
                this->dAplus_dqx[gp] * (q - q_hat) - this->dAminus_dqx[gp] * (q_inf - q_hat);
            column(dB_dq_hat, SWE::Variables::qy) +=
                this->dAplus_dqy[gp] * (q - q_hat) - this->dAminus_dqy[gp] * (q_inf - q_hat);

            column(edge_internal.delta_hat_global_kernel_at_gp, gp) = flatten<double>(dB_dq_hat);
            column(edge_internal.rhs_global_kernel_at_gp, gp) =
                this->Aplus[gp] * (q - q_hat) - this->Aminus[gp] * (q_inf - q_hat);
        }

        for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
            for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
                submatrix(edge_internal.delta_hat_global,
                          SWE::n_variables * dof_i,
                          SWE::n_variables * dof_j,
                          SWE::n_variables,
                          SWE::n_variables) =
                    reshape<double, SWE::n_variables>(
                        edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
            }

            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
                -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }

        solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

        edge_state.q_hat += reshape<double, SWE::n_variables, SO::ColumnMajor>(edge_internal.rhs_global,
                                                                               edge_bound.edge_data.get_ndof());

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-12) {
            break;
        }
    }
}

template <typename StepperType, typename EdgeBoundaryType>
void Function::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    get_Aplus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aplus);
    get_dAplus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dze);
    get_dAplus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqx);
    get_dAplus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqy);

    get_Aminus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aminus);
    get_dAminus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dze);
    get_dAminus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqx);
    get_dAminus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqy);

    double t = stepper.GetTimeAtCurrentStage() + stepper.GetDT();

    HybMatrix<double, SWE::n_variables> q_ex = edge_bound.boundary.ComputeFgp([t](Point<2>& pt) {
        double ze = 0.0;
        double qx = 0.0;
        double qy = 0.0;

        Utilities::ignore(ze, qx, qy);

        if (t <= 3.0) {
            ze = cos(PI * t) - 1.0;
        } else {
            ze = -2.0;
        }

        // StatVector<double, SWE::n_variables> q{ze, qx, qy};
        StatVector<double, SWE::n_variables> q(SWE::ic_q(t, pt));

        return q;
    });

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        auto q     = column(boundary.q_at_gp, gp);
        auto q_hat = column(q_hat_at_gp, gp);
        auto q_inf = column(q_ex, gp);

        StatMatrix<double, SWE::n_variables, SWE::n_variables> dB_dq_hat = this->Aminus[gp] - this->Aplus[gp];

        column(dB_dq_hat, SWE::Variables::ze) +=
            this->dAplus_dze[gp] * (q - q_hat) - this->dAminus_dze[gp] * (q_inf - q_hat);
        column(dB_dq_hat, SWE::Variables::qx) +=
            this->dAplus_dqx[gp] * (q - q_hat) - this->dAminus_dqx[gp] * (q_inf - q_hat);
        column(dB_dq_hat, SWE::Variables::qy) +=
            this->dAplus_dqy[gp] * (q - q_hat) - this->dAminus_dqy[gp] * (q_inf - q_hat);

        column(boundary.delta_global_kernel_at_gp, gp)          = flatten<double>(this->Aplus[gp]);
        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = flatten<double>(dB_dq_hat);
        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            this->Aplus[gp] * (q - q_hat) - this->Aminus[gp] * (q_inf - q_hat);
    }
}
}
}
}

#endif
