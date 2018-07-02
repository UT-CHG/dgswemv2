#ifndef RKDG_SWE_DATA_SLOPE_LIMIT_HPP
#define RKDG_SWE_DATA_SLOPE_LIMIT_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nbound)
        : surface_normal(nbound, std::vector<double>(2)),
          midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          alpha_1(nbound),
          alpha_2(nbound),
          r_sq(nbound),
          ze_lin(nbound),
          qx_lin(nbound),
          qy_lin(nbound),
          ze_at_vrtx(nbound),
          qx_at_vrtx(nbound),
          qy_at_vrtx(nbound),
          ze_at_midpts(nbound),
          qx_at_midpts(nbound),
          qy_at_midpts(nbound),
          bath_at_midpts(nbound),
          wet_neigh(nbound),
          ze_at_baryctr_neigh(nbound),
          qx_at_baryctr_neigh(nbound),
          qy_at_baryctr_neigh(nbound) {}

    Array2D<double> surface_normal;

    Point<2> baryctr_coord;
    std::vector<Point<2>> midpts_coord;
    std::vector<Point<2>> baryctr_coord_neigh;

    std::vector<double> alpha_1;
    std::vector<double> alpha_2;
    std::vector<double> r_sq;

    std::vector<double> ze_lin;
    std::vector<double> qx_lin;
    std::vector<double> qy_lin;

    double ze_at_baryctr;
    double qx_at_baryctr;
    double qy_at_baryctr;
    double bath_at_baryctr;

    std::vector<double> ze_at_vrtx;
    std::vector<double> qx_at_vrtx;
    std::vector<double> qy_at_vrtx;

    std::vector<double> ze_at_midpts;
    std::vector<double> qx_at_midpts;
    std::vector<double> qy_at_midpts;
    std::vector<double> bath_at_midpts;

    std::vector<bool> wet_neigh;
    std::vector<double> ze_at_baryctr_neigh;
    std::vector<double> qx_at_baryctr_neigh;
    std::vector<double> qy_at_baryctr_neigh;

    std::vector<double> w_midpt_char = std::vector<double>(3);
    Array2D<double> w_baryctr_char   = Array2D<double>(3, std::vector<double>(3));

    std::vector<double> delta_char = std::vector<double>(3);
    Array2D<double> delta          = Array2D<double>(3, std::vector<double>(3));

    Array2D<double> L = Array2D<double>(3, std::vector<double>(3));
    Array2D<double> R = Array2D<double>(3, std::vector<double>(3));

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void SlopeLimit::serialize(Archive& ar, unsigned) {
    // clang-format off
    ar  & surface_normal
        & baryctr_coord
        & midpts_coord
        & baryctr_coord_neigh
        & alpha_1
        & alpha_2
        & r_sq
        & ze_lin
        & qx_lin
        & qy_lin
        & ze_at_baryctr
        & qx_at_baryctr
        & qy_at_baryctr
        & bath_at_baryctr
        & ze_at_vrtx
        & qx_at_vrtx
        & qy_at_vrtx
        & ze_at_midpts
        & qx_at_midpts
        & qy_at_midpts
        & bath_at_midpts
        & wet_neigh
        & ze_at_baryctr_neigh
        & qx_at_baryctr_neigh
        & qy_at_baryctr_neigh
        & w_midpt_char
        & w_baryctr_char
        & delta_char
        & delta
        & L
        & R;
    // clang-format on
}
#endif
}
}

#endif
