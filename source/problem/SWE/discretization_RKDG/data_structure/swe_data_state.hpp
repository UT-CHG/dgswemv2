#ifndef SWE_DATA_STATE_HPP
#define SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct State {
    State() = default;
    State(const uint ndof)
        : ze(ndof),
          qx(ndof),
          qy(ndof),
          bath(ndof),
          rhs_ze(ndof),
          rhs_qx(ndof),
          rhs_qy(ndof),
          solution_ze(ndof),
          solution_qx(ndof),
          solution_qy(ndof) {}

    std::vector<double> ze;
    std::vector<double> qx;
    std::vector<double> qy;
    std::vector<double> bath;

    std::vector<double> rhs_ze;
    std::vector<double> rhs_qx;
    std::vector<double> rhs_qy;

    std::vector<double> solution_ze;
    std::vector<double> solution_qx;
    std::vector<double> solution_qy;
};
}
}

#endif
