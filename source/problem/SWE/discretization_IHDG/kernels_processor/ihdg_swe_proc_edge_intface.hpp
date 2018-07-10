#ifndef IHDG_SWE_PROC_EDGE_INTFACE_HPP
#define IHDG_SWE_PROC_EDGE_INTFACE_HPP

#include <eigen3/Eigen/Dense>

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
void Problem::prepare_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_intface_LF(edge_int);

    // Add dtau/dq * del_q - tau * del_q term to dF_hat_dq_hat
    // and tau * del_q term to dF_hat_dq
    add_dF_hat_tau_terms_intface_LF(edge_int);

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);
}
}
}

#endif