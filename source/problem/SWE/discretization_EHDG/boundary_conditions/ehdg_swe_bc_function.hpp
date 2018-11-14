#ifndef EHDG_SWE_BC_FUNCTION_HPP
#define EHDG_SWE_BC_FUNCTION_HPP

namespace SWE {
namespace EHDG {
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
