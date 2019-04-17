#ifndef EHDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP
#define EHDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP

#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
class DistributedLevee {
  public:
    DBDataExchanger exchanger;

  private:
    double H_tolerance = 0.01;

    DynRowVector<double> H_barrier;
    DynRowVector<double> C_subcritical;
    DynRowVector<double> C_supercritical;

    DynRowVector<double> H_bar_gp;
    DynRowVector<double> C_subcrit_gp;
    DynRowVector<double> C_supercrit_gp;

    BC::Land land_boundary;

  public:
    DistributedLevee() = default;
    DistributedLevee(const DBDataExchanger& exchanger, const std::vector<LeveeInput>& levee_input);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(EdgeDistributedType& edge_dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);
};

DistributedLevee::DistributedLevee(const DBDataExchanger& exchanger, const std::vector<LeveeInput>& levee_input)
    : exchanger(exchanger) {
    uint n_nodes = levee_input.size();

    this->H_barrier.resize(n_nodes);
    this->C_subcritical.resize(n_nodes);
    this->C_supercritical.resize(n_nodes);

    for (uint node = 0; node < n_nodes; ++node) {
        this->H_barrier[node]       = levee_input[node].H_barrier;
        this->C_subcritical[node]   = levee_input[node].C_subcritical;
        this->C_supercritical[node] = levee_input[node].C_supercritical;
    }
}

template <typename DistributedBoundaryType>
void DistributedLevee::Initialize(DistributedBoundaryType& dbound) {
    this->H_bar_gp       = dbound.ComputeBoundaryNodalUgp(this->H_barrier);
    this->C_subcrit_gp   = dbound.ComputeBoundaryNodalUgp(this->C_subcritical);
    this->C_supercrit_gp = dbound.ComputeBoundaryNodalUgp(this->C_supercritical);
}

template <typename EdgeDistributedType>
void DistributedLevee::ComputeGlobalKernels(EdgeDistributedType& edge_dbound) {
    // Something to implement in the future
}

template <typename EdgeDistributedType>
void DistributedLevee::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    // Something to implement in the future
}
}
}
}

#endif