#ifndef EHDG_SWE_IS_LEVEE_HPP
#define EHDG_SWE_IS_LEVEE_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_EHDG/stabilization_parameters/ehdg_swe_stabilization_params.hpp"

namespace SWE {
namespace EHDG {
namespace ISP {
class Levee {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFlux(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Levee::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    // Something to implement in the future
}

template <typename EdgeInterfaceType>
void Levee::ComputeNumericalFlux(EdgeInterfaceType& edge_int) {
    // Something to implement in the future
}
}
}
}

#endif