#ifndef SWE_DATA_STATE_HPP
#define SWE_DATA_STATE_HPP

#include "../../../general_definitions.hpp"

namespace SWE {
struct State {
    State() = default;
    State(const uint ndof) : ze(ndof), qx(ndof), qy(ndof), bath(ndof), rhs_ze(ndof), rhs_qx(ndof), rhs_qy(ndof) {}

    std::vector<double> ze;
    std::vector<double> qx;
    std::vector<double> qy;
    std::vector<double> bath;

    std::vector<double> rhs_ze;
    std::vector<double> rhs_qx;
    std::vector<double> rhs_qy;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void State::serialize(Archive& ar, unsigned) {
    ar & ze & qx & qy & bath & rhs_ze & rhs_qx & rhs_qy;
}
#endif
}

#endif
