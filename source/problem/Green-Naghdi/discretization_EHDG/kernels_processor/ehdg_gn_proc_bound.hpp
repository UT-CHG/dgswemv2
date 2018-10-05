#ifndef EHDG_GN_PROC_BOUND_HPP
#define EHDG_GN_PROC_BOUND_HPP

namespace GN {
namespace EHDG {
template <typename BoundaryType>
void Problem::global_swe_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    SWE::EHDG::Problem::global_boundary_kernel(stepper, bound);
}

template <typename BoundaryType>
void Problem::local_swe_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    SWE::EHDG::Problem::local_boundary_kernel(stepper, bound);
}
}
}

#endif
