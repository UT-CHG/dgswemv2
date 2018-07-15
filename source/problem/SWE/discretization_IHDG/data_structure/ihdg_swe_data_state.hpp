#ifndef IHDG_SWE_DATA_STATE_HPP
#define IHDG_SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct State {
    State() = default;
    State(const uint ndof) : q(ndof), bath(ndof), rhs(ndof), solution(ndof) {}

    std::vector<StatVector<double, SWE::n_variables>> q;
    std::vector<double> bath;
    std::vector<StatVector<double, SWE::n_variables>> rhs;
    std::vector<StatVector<double, SWE::n_variables>> solution;
};
}
}

#endif
