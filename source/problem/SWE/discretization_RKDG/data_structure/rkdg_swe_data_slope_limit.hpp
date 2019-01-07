#ifndef RKDG_SWE_DATA_SLOPE_LIMIT_HPP
#define RKDG_SWE_DATA_SLOPE_LIMIT_HPP

namespace SWE {
namespace RKDG {
struct SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nvrtx, const uint nbound)
        : midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          median(nbound),
          alpha(nbound),
          r_sq(nbound),
          q_lin(SWE::n_variables, nvrtx),
          q_at_vrtx(SWE::n_variables, nvrtx),
          q_at_midpts(SWE::n_variables, nbound),
          bath_at_vrtx(nvrtx),
          bath_at_midpts(nbound),
          wet_neigh(nbound),
          q_at_baryctr_neigh(nbound),
          delta(SWE::n_variables, nbound) {}

    Point<2> baryctr_coord;
    std::vector<Point<2>> midpts_coord;
    std::vector<Point<2>> baryctr_coord_neigh;

    AlignedVector<StatVector<double, SWE::n_dimensions>> median;
    AlignedVector<StatVector<double, SWE::n_dimensions>> alpha;
    std::vector<double> r_sq;

    HybMatrix<double, SWE::n_variables> q_lin;

    StatVector<double, SWE::n_variables> q_at_baryctr;
    HybMatrix<double, SWE::n_variables> q_at_vrtx;
    HybMatrix<double, SWE::n_variables> q_at_midpts;

    double bath_at_baryctr;
    DynRowVector<double> bath_at_vrtx;
    DynRowVector<double> bath_at_midpts;

    std::vector<bool> wet_neigh;
    AlignedVector<StatVector<double, SWE::n_variables>> q_at_baryctr_neigh;

    HybMatrix<double, SWE::n_variables> delta;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & surface_normal
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
            & delta;
        // clang-format on
    }
#endif
};
}
}

#endif
