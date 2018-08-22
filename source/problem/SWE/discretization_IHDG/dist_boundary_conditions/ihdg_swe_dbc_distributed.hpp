#ifndef IHDG_SWE_DBC_DISTRIBUTED_HPP
#define IHDG_SWE_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace IHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound) {} /*nothing to initialize*/

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {}
}
}
}

#endif