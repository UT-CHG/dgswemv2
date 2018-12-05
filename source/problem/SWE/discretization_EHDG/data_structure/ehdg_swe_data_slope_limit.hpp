#ifndef EHDG_SWE_DATA_SLOPE_LIMIT_HPP
#define EHDG_SWE_DATA_SLOPE_LIMIT_HPP

namespace SWE {
namespace EHDG {
struct SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nvrtx, const uint nbound)
        : surface_normal(nbound),
          midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          alpha_1(nbound),
          alpha_2(nbound),
          r_sq(nbound),
          q_lin(SWE::n_variables, nvrtx),
          q_at_vrtx(SWE::n_variables, nvrtx),
          q_at_midpts(SWE::n_variables, nbound),
          bath_at_vrtx(nvrtx),
          bath_at_midpts(nbound),
          wet_neigh(nbound),
          q_at_baryctr_neigh(nbound) {}

    AlignedVector<StatVector<double, SWE::n_dimensions>> surface_normal;

    Point<2> baryctr_coord;
    std::vector<Point<2>> midpts_coord;
    std::vector<Point<2>> baryctr_coord_neigh;

    std::vector<double> alpha_1;
    std::vector<double> alpha_2;
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

    StatMatrix<double, SWE::n_variables, SWE::n_variables> delta_char;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> delta;
};
}
}

#endif
