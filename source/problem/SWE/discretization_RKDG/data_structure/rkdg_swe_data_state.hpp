#ifndef RKDG_SWE_DATA_STATE_HPP
#define RKDG_SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct State {
    State() = default;
    State(const uint ndof) : q(ndof), bath(ndof), rhs(ndof), solution(ndof) {}

    std::vector<StatVector<double, SWE::n_variables>> q;
    std::vector<double> bath;
    std::vector<StatVector<double, SWE::n_variables>> rhs;
    std::vector<StatVector<double, SWE::n_variables>> solution;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        /*ar  & q
            & bath
            & rhs
            & solution;*/
        // clang-format on
    }
#endif
};
}
}

#endif
