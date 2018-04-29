#ifndef SWE_INTFACE_SPEC_LEVEE_HPP
#define SWE_INTFACE_SPEC_LEVEE_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace IS {
class Levee {
  public:
    template <typename InterfaceType>
    void ComputeFlux(const Stepper& stepper, InterfaceType& intface);
};

template <typename InterfaceType>
void Levee::ComputeFlux(const Stepper& stepper, InterfaceType& intface) {}
}
}

#endif