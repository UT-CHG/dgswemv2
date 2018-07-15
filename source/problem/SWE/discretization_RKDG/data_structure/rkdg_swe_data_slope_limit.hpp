#ifndef RKDG_SWE_DATA_SLOPE_LIMIT_HPP
#define RKDG_SWE_DATA_SLOPE_LIMIT_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nvrtx, const uint nbound)
        : surface_normal(nbound, std::vector<double>(2)),
          midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          alpha_1(nbound),
          alpha_2(nbound),
          r_sq(nbound),
          q_lin(nvrtx),
          q_at_vrtx(nvrtx),
          q_at_midpts(nbound),
          bath_at_vrtx(nvrtx),
          bath_at_midpts(nbound),
          wet_neigh(nbound),
          q_at_baryctr_neigh(nbound) {}

    Array2D<double> surface_normal;

    Point<2> baryctr_coord;
    std::vector<Point<2>> midpts_coord;
    std::vector<Point<2>> baryctr_coord_neigh;

    std::vector<double> alpha_1;
    std::vector<double> alpha_2;
    std::vector<double> r_sq;

    std::vector<StatVector<double, SWE::n_variables>> q_lin;

    StatVector<double, SWE::n_variables> q_at_baryctr;
    std::vector<StatVector<double, SWE::n_variables>> q_at_vrtx;
    std::vector<StatVector<double, SWE::n_variables>> q_at_midpts;

    double bath_at_baryctr;
    std::vector<double> bath_at_vrtx;
    std::vector<double> bath_at_midpts;

    std::vector<bool> wet_neigh;
    std::vector<StatVector<double, SWE::n_variables>> q_at_baryctr_neigh;

    StatMatrix<double, SWE::n_variables, SWE::n_variables> L;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> R;

    StatVector<double, SWE::n_variables> w_midpt_char;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> w_baryctr_char;

    StatVector<double, SWE::n_variables> delta_char;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> delta;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        /*ar  & surface_normal
            & baryctr_coord
            & midpts_coord
            & baryctr_coord_neigh
            & alpha_1
            & alpha_2
            & r_sq
            & q_lin
            & q_at_baryctr
            & q_at_vrtx
            & q_at_midpts
            & bath_at_baryctr
            & bath_at_vrtx
            & bath_at_midpts
            & wet_neigh
            & q_at_baryctr_neigh
            & w_midpt_char
            & w_baryctr_char
            & delta_char
            & delta
            & L
            & R;*/
        // clang-format on
    }
#endif
};
}
}

#endif
