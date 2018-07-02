#ifndef RKDG_SWE_DATA_STATE_HPP
#define RKDG_SWE_DATA_STATE_HPP

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

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void State::serialize(Archive& ar, unsigned) {
    // clang-format off
    ar  & ze
        & qx
        & qy
        & bath
        & rhs_ze
        & rhs_qx
        & rhs_qy
        & solution_ze
        & solution_qx
        & solution_qy;
    // clang-format on
}
#endif
}
}

#endif
