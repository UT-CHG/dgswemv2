#ifndef IHDG_SWE_BC_TIDE_HPP
#define IHDG_SWE_BC_TIDE_HPP

#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"

namespace SWE {
namespace IHDG {
namespace BC {
class Tide {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;

    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<DynRowVector<double>> amplitude;
    std::vector<DynRowVector<double>> phase;

    std::vector<DynRowVector<double>> amplitude_gp;
    std::vector<DynRowVector<double>> phase_gp;

    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> Aplus;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dqy;

    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> Aminus;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dqy;

  public:
    Tide() = default;
    Tide(const std::vector<TideNode>& tide_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

Tide::Tide(const std::vector<TideNode>& tide_input) {
    this->frequency    = tide_input[0].frequency;
    this->forcing_fact = tide_input[0].forcing_fact;
    this->equilib_arg  = tide_input[0].equilib_arg;

    uint n_contituents = this->frequency.size();
    uint n_nodes       = tide_input.size();

    this->amplitude.resize(n_contituents);
    this->phase.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude[con].resize(n_nodes);
        this->phase[con].resize(n_nodes);

        for (uint node = 0; node < n_nodes; ++node) {
            this->amplitude[con][node] = tide_input[node].amplitude[con];
            this->phase[con][node]     = tide_input[node].phase[con];
        }
    }
}

template <typename BoundaryType>
void Tide::Initialize(BoundaryType& bound) {
    uint ngp           = bound.data.get_ngp_boundary(bound.bound_id);
    uint n_contituents = this->frequency.size();

    this->q_ex.resize(SWE::n_variables, ngp);

    this->amplitude_gp.resize(n_contituents);
    this->phase_gp.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude_gp[con] = bound.ComputeBoundaryNodalUgp(this->amplitude[con]);
        this->phase_gp[con]     = bound.ComputeBoundaryNodalUgp(this->phase[con]);
    }

    this->Aplus.resize(ngp);
    this->dAplus_dze.resize(ngp);
    this->dAplus_dqx.resize(ngp);
    this->dAplus_dqy.resize(ngp);

    this->Aminus.resize(ngp);
    this->dAminus_dze.resize(ngp);
    this->dAminus_dqx.resize(ngp);
    this->dAminus_dqy.resize(ngp);
}

template <typename StepperType, typename EdgeBoundaryType>
void Tide::ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& bound = edge_bound.boundary;

    auto& state    = bound.data.state[0];
    auto& boundary = bound.data.boundary[bound.bound_id];

    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    boundary.q_at_gp = bound.ComputeUgp(state.q);

    set_constant(edge_state.q_hat, 0.0);

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        get_Aplus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aplus);
        get_dAplus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dze);
        get_dAplus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqx);
        get_dAplus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqy);

        get_Aminus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aminus);
        get_dAminus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dze);
        get_dAminus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqx);
        get_dAminus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqy);

        set_constant(this->q_ex, 0.0);

        for (uint con = 0; con < this->frequency.size(); ++con) {
            for (uint gp = 0; gp < columns(boundary.q_at_gp); ++gp) {
                this->q_ex(SWE::Variables::ze, gp) += stepper.GetRamp() * this->forcing_fact[con] *
                                                      this->amplitude_gp[con][gp] *
                                                      cos(this->frequency[con] * stepper.GetTimeAtCurrentStage() +
                                                          this->equilib_arg[con] - this->phase_gp[con][gp]);
            }
        }

        row(this->q_ex, SWE::Variables::qx) = row(boundary.q_at_gp, SWE::Variables::qx);
        row(this->q_ex, SWE::Variables::qy) = row(boundary.q_at_gp, SWE::Variables::qy);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            auto q     = column(boundary.q_at_gp, gp);
            auto q_hat = column(edge_internal.q_hat_at_gp, gp);
            auto q_inf = column(this->q_ex, gp);

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
void Tide::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
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

    set_constant(this->q_ex, 0.0);

    for (uint con = 0; con < this->frequency.size(); ++con) {
        for (uint gp = 0; gp < columns(boundary.q_at_gp); ++gp) {
            this->q_ex(SWE::Variables::ze, gp) += stepper.GetRampNext() * this->forcing_fact[con] *
                                                  this->amplitude_gp[con][gp] *
                                                  cos(this->frequency[con] * stepper.GetTimeAtNextStage() +
                                                      this->equilib_arg[con] - this->phase_gp[con][gp]);
        }
    }

    row(this->q_ex, SWE::Variables::qx) = row(boundary.q_at_gp, SWE::Variables::qx);
    row(this->q_ex, SWE::Variables::qy) = row(boundary.q_at_gp, SWE::Variables::qy);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        auto q     = column(boundary.q_at_gp, gp);
        auto q_hat = column(q_hat_at_gp, gp);
        auto q_inf = column(this->q_ex, gp);

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
