#ifndef EHDG_GN_PROC_EDGE_DBOUND_HPP
#define EHDG_GN_PROC_EDGE_DBOUND_HPP

namespace GN {
namespace EHDG {
template <typename EdgeDistributedType>
void Problem::global_swe_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    SWE::EHDG::Problem::global_edge_distributed_kernel(stepper, edge_dbound);
}
}
}

#endif