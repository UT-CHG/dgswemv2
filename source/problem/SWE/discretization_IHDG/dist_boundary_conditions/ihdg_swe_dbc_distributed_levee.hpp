#ifndef IHDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP
#define IHDG_SWE_DBC_DISTRIBUTED_LEVEE_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace IHDG {
namespace DBC {
class DistributedLevee {
  public:
    DBDataExchanger exchanger;

  public:
    DistributedLevee() = default;
    DistributedLevee(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound) {} /*nothing to initialize*/

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound);
};

DistributedLevee::DistributedLevee(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename EdgeDistributedType>
void DistributedLevee::ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    // Something to implement in the future
}
}
}
}

#endif