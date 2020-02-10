#ifndef RKDG_SWE_BC_TIDE_HPP
#define RKDG_SWE_BC_TIDE_HPP

namespace SWE {
namespace RKDG {
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

  public:
    Tide() = default;
    Tide(const std::vector<TideNode>& tide_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename BoundaryType>
    void ComputeFlux(const StepperType& stepper, BoundaryType& bound);

    template <typename StepperType, typename BoundaryType>
    void ComputeBedFlux(const StepperType& stepper, BoundaryType& bound);
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
}

template <typename StepperType, typename BoundaryType>
void Tide::ComputeFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];

    set_constant(this->q_ex, 0.0);

    for (uint con = 0; con < this->frequency.size(); ++con) {
        for (uint gp = 0; gp < columns(boundary.q_at_gp); ++gp) {
            this->q_ex(SWE::Variables::ze, gp) += stepper.GetRamp() * this->forcing_fact[con] *
                                                  this->amplitude_gp[con][gp] *
                                                  cos(this->frequency[con] * stepper.GetTimeAtCurrentStage() +
                                                      (this->equilib_arg[con] - this->phase_gp[con][gp]) * PI / 180);
        }
    }

    row(this->q_ex, SWE::Variables::qx) = row(boundary.q_at_gp, SWE::Variables::qx);
    row(this->q_ex, SWE::Variables::qy) = row(boundary.q_at_gp, SWE::Variables::qy);
    set_constant(row(this->q_ex, SWE::Variables::hc), 0.0);  // TODO

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        HLL_flux(Global::g,
                 column(boundary.q_at_gp, gp),
                 column(this->q_ex, gp),
                 column(boundary.aux_at_gp, gp),
                 column(boundary.aux_at_gp, gp),
                 column(bound.surface_normal, gp),
                 column(boundary.F_hat_at_gp, gp));
    }
}

template <typename StepperType, typename BoundaryType>
void Tide::ComputeBedFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        const double un = roe_un(column(boundary.q_at_gp, gp),
                                 column(this->q_ex, gp),
                                 column(boundary.aux_at_gp, gp),
                                 column(boundary.aux_at_gp, gp),
                                 column(bound.surface_normal, gp));
        if (Utilities::almost_equal(un, 0.0)) {
            boundary.qb_hat_at_gp[gp] = 0.0;
        } else if (un > 0.0) {
            boundary.qb_hat_at_gp[gp] = transpose(column(bound.surface_normal, gp)) *
                                        bed_flux(column(boundary.q_at_gp, gp), column(boundary.aux_at_gp, gp));
        } else if (un < 0.0) {
            boundary.qb_hat_at_gp[gp] = transpose(column(bound.surface_normal, gp)) *
                                        bed_flux(column(this->q_ex, gp), column(boundary.aux_at_gp, gp));
        }
    }
}
}
}
}

#endif
