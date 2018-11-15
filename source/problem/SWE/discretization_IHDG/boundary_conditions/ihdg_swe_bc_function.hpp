#ifndef IHDG_SWE_BC_FUNCTION_HPP
#define IHDG_SWE_BC_FUNCTION_HPP

namespace SWE {
namespace IHDG {
namespace BC {
class Function {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {}

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {}

    template <typename EdgeBoundaryType>
    void ComputeNumericalFlux(EdgeBoundaryType& edge_bound) {}
};
}
}
}

#endif