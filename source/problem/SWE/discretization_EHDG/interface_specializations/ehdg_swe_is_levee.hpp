#ifndef EHDG_SWE_IS_LEVEE_HPP
#define EHDG_SWE_IS_LEVEE_HPP

namespace SWE {
namespace EHDG {
namespace ISP {
class Levee {
  private:
    // double H_tolerance = 0.01;

    DynRowVector<double> H_barrier;
    DynRowVector<double> C_subcritical;
    DynRowVector<double> C_supercritical;

    DynRowVector<double> H_bar_gp;
    DynRowVector<double> C_subcrit_gp;
    DynRowVector<double> C_supercrit_gp;

    // BC::Land land_boundary;

  public:
    Levee() = default;
    Levee(const std::vector<LeveeInput>& levee_input);

    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFlux(EdgeInterfaceType& edge_int);

    template <typename InterfaceType>
    void ComputeBedFlux(InterfaceType& intface) {
        abort();  // no implementation
    }
};

Levee::Levee(const std::vector<LeveeInput>& levee_input) {
    uint n_nodes = levee_input.size();

    this->H_barrier.resize(n_nodes);
    this->C_subcritical.resize(n_nodes);
    this->C_supercritical.resize(n_nodes);

    for (uint node = 0; node < n_nodes; ++node) {
        this->H_barrier[node]       = levee_input[node].H_barrier;
        this->C_subcritical[node]   = levee_input[node].C_subcritical;
        this->C_supercritical[node] = levee_input[node].C_supercritical;
    }
}

template <typename InterfaceType>
void Levee::Initialize(InterfaceType& intface) {
    this->H_bar_gp       = intface.ComputeBoundaryNodalUgpIN(this->H_barrier);
    this->C_subcrit_gp   = intface.ComputeBoundaryNodalUgpIN(this->C_subcritical);
    this->C_supercrit_gp = intface.ComputeBoundaryNodalUgpIN(this->C_supercritical);
}

template <typename EdgeInterfaceType>
void Levee::ComputeNumericalFlux(EdgeInterfaceType& edge_int) {
    abort();  // no implementation
}
}
}
}

#endif