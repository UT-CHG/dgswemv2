#ifndef EHDG_GN_BC_FLOW_HPP
#define EHDG_GN_BC_FLOW_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"

namespace GN {
namespace EHDG {
namespace BC {
class Flow {
  private:
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    Array2D<double> amplitude;
    Array2D<double> phase;

    Array2D<double> amplitude_gp;
    Array2D<double> phase_gp;

  public:
    Flow() = default;
    Flow(const std::vector<FlowInput>& flow_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    /* something to be implemented in future */

    template <typename EdgeBoundaryType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {}

    template <typename EdgeBoundaryType>
    void ComputeNumericalFlux(EdgeBoundaryType& edge_bound) {}

    template <typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {}
};

Flow::Flow(const std::vector<FlowInput>& flow_input) {
    this->frequency    = flow_input[0].frequency;
    this->forcing_fact = flow_input[0].forcing_fact;
    this->equilib_arg  = flow_input[0].equilib_arg;

    uint n_contituents = this->frequency.size();
    uint n_nodes       = flow_input.size();

    this->amplitude.resize(n_contituents);
    this->phase.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude[con].resize(n_nodes);
        this->phase[con].resize(n_nodes);

        for (uint node = 0; node < n_nodes; ++node) {
            this->amplitude[con][node] = flow_input[node].amplitude[con];
            this->phase[con][node]     = flow_input[node].phase[con];
        }
    }
}

template <typename BoundaryType>
void Flow::Initialize(BoundaryType& bound) {
    uint n_contituents = this->frequency.size();

    this->amplitude_gp.resize(n_contituents);
    this->phase_gp.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude_gp[con].resize(bound.data.get_ngp_boundary(bound.bound_id));
        this->phase_gp[con].resize(bound.data.get_ngp_boundary(bound.bound_id));

        // bound.ComputeBoundaryNodalUgp(this->amplitude[con], this->amplitude_gp[con]);
        // bound.ComputeBoundaryNodalUgp(this->phase[con], this->phase_gp[con]);
    }
}
}
}
}

#endif