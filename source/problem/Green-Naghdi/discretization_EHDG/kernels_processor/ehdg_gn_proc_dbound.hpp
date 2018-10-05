#ifndef EHDG_GN_PROC_DBOUND_HPP
#define EHDG_GN_PROC_DBOUND_HPP

namespace GN {
namespace EHDG {
template <typename DistributedBoundaryType>
void Problem::global_swe_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    SWE::EHDG::Problem::global_distributed_boundary_kernel(stepper, dbound);
}

template <typename DistributedBoundaryType>
void Problem::local_swe_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    SWE::EHDG::Problem::local_distributed_boundary_kernel(stepper, dbound);
}
}
}

#endif
