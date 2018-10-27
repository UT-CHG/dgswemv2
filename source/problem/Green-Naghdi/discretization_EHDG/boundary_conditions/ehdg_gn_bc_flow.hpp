#ifndef EHDG_GN_BC_FLOW_HPP
#define EHDG_GN_BC_FLOW_HPP

namespace GN {
namespace EHDG {
namespace BC {
class Flow : public SWE::EHDG::BC::Flow {
  public:
    Flow() = default;
    Flow(const std::vector<SWE::FlowInput>& flow_input) : SWE::EHDG::BC::Flow(flow_input) {}

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound) {}
};
}
}
}

#endif