#ifndef EHDG_GN_PROC_INTFACE_HPP
#define EHDG_GN_PROC_INTFACE_HPP

namespace GN {
namespace EHDG {
template <typename InterfaceType>
void Problem::global_swe_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    SWE::EHDG::Problem::global_interface_kernel(stepper, intface);
}

template <typename InterfaceType>
void Problem::local_swe_interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    SWE::EHDG::Problem::local_interface_kernel(stepper, intface);
}
}
}

#endif
