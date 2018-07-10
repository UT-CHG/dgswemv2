#ifndef IHDG_SWE_PROC_EDGE_INTFACE_HPP
#define IHDG_SWE_PROC_EDGE_INTFACE_HPP

#include <eigen3/Eigen/Dense>

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
void Problem::prepare_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    edge_int.interface.specialization.ComputeNumericalFlux(edge_int);

    add_dF_hat_tau_terms_intface_LF(edge_int);

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);
}
}
}

#endif