#ifndef EHDG_GN_IS_LEVEE_HPP
#define EHDG_GN_IS_LEVEE_HPP

namespace GN {
namespace EHDG {
namespace ISP {
class Levee : public SWE_SIM::ISP::Levee {
  public:
    Levee() = default;
    Levee(const std::vector<SWE::LeveeInput>& levee_input) : SWE_SIM::ISP::Levee(levee_input){};

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