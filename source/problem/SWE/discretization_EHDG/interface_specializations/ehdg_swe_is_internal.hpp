#ifndef EHDG_SWE_IS_INTERNAL_HPP
#define EHDG_SWE_IS_INTERNAL_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_EHDG/numerical_fluxes/ehdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace EHDG {
namespace IS {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/
};
}
}
}

#endif