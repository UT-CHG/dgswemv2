#ifndef EHDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP
#define EHDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
class DistributedLevee {
  public:
    DBDataExchanger exchanger;

  public:
    DistributedLevee() = default;
    DistributedLevee(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);
};

DistributedLevee::DistributedLevee(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void DistributedLevee::Initialize(DistributedBoundaryType& dbound) {
    // Something to implement in the future
}

template <typename EdgeDistributedType>
void DistributedLevee::ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
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