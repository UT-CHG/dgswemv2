#ifndef EHDG_GN_BC_TIDE_HPP
#define EHDG_GN_BC_TIDE_HPP

namespace GN {
namespace EHDG {
namespace BC {
class Tide : public SWE::EHDG::BC::Tide {
  public:
    Tide() = default;
    Tide(const std::vector<SWE::TideInput>& tide_input) : SWE::EHDG::BC::Tide(tide_input) {}

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound) {}
};
}
}
}

#endif
