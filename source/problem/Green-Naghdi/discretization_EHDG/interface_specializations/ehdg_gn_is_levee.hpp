#ifndef EHDG_GN_IS_LEVEE_HPP
#define EHDG_GN_IS_LEVEE_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/stabilization_parameters/ehdg_gn_stabilization_params.hpp"

namespace GN {
namespace EHDG {
namespace ISP {
class Levee : public SWE::EHDG::ISP::Levee {
  public:
    Levee() = default;
    Levee(const std::vector<SWE::LeveeInput>& levee_input) : SWE::EHDG::ISP::Levee(levee_input){};

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernelsDC(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Levee::ComputeGlobalKernelsDC(EdgeInterfaceType& edge_int) {
    // Something to implement in the future
}
}
}
}

#endif