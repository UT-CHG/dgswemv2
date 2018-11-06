#ifndef EHDG_GN_DBC_DISTRIBUTED_LEVEE_HPP
#define EHDG_GN_DBC_DISTRIBUTED_LEVEE_HPP

#include "communication/db_data_exchanger.hpp"

namespace GN {
namespace EHDG {
namespace DBC {
class DistributedLevee : public SWE_SIM::DBC::DistributedLevee {
  public:
    DistributedLevee() = default;
    DistributedLevee(const DBDataExchanger& exchanger, const std::vector<SWE::LeveeInput>& levee_input)
        : SWE_SIM::DBC::DistributedLevee(exchanger, levee_input) {}

    template <typename EdgeDistributedType>
    void ComputeGlobalKernelsDC(EdgeDistributedType& edge_dbound);
};

template <typename EdgeDistributedType>
void DistributedLevee::ComputeGlobalKernelsDC(EdgeDistributedType& edge_dbound) {
    // Something to implement in the future
}
}
}
}

#endif