#ifndef SWE_DATA_SLOPE_LIMIT_HPP
#define SWE_DATA_SLOPE_LIMIT_HPP

namespace SWE {
struct SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nvrtx, const uint nbound)
        : hdif_at_gp(nbound),
          midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          median(nbound),
          alpha(nbound),
          a_elem(2 * nbound),
          r_sq(nbound),
          A_inv(nbound),
          inv_theta_r(nbound),
          d_r(nbound),
          alpha_r(nbound),
          w_r(nbound),
          q_lin(SWE::n_variables, nvrtx),
          bath_lin(nvrtx),
          q_at_vrtx(SWE::n_variables, nvrtx),
          q_at_midpts(SWE::n_variables, nbound),
          bath_at_vrtx(nvrtx),
          bath_at_midpts(nbound),
          wet_neigh(nbound),
          q_at_baryctr_neigh(nbound),
          bath_at_baryctr_neigh(nbound),
          delta(SWE::n_variables, nbound),
          bath_delta(nbound),
          dq_r(nbound),
          dbath_r(SWE::n_dimensions, nbound) {
        // *** //
        this->T = DynMatrix<double>(nvrtx, nvrtx);  // This is for a triangular element!
        set_constant(this->T, 1.0);
        for (uint vrtx = 0; vrtx < nvrtx; ++vrtx) {
            this->T(vrtx, vrtx) = -1.0;
        }
    }

    std::vector<double> lengths;
    double radius;
    bool troubled    = false;
    double perimeter = 0.0;
    double I         = 0.0;
    std::vector<DynRowVector<double>> hdif_at_gp;

    DynMatrix<double> T;

    Point<2> baryctr_coord;
    AlignedVector<Point<2>> midpts_coord;
    AlignedVector<Point<2>> baryctr_coord_neigh;

    // COCKBURN-SHU
    AlignedVector<StatVector<double, SWE::n_dimensions>> median;
    AlignedVector<StatVector<double, SWE::n_dimensions>> alpha;
    std::vector<uint> a_elem;
    std::vector<double> r_sq;

    // XU ETAL
    AlignedVector<StatMatrix<double, SWE::n_dimensions, SWE::n_dimensions>> A_inv;
    std::vector<double> inv_theta_r;
    std::vector<double> d_r;
    std::vector<double> alpha_r;
    DynVector<double> w_r;

    HybMatrix<double, SWE::n_variables> q_lin;
    DynRowVector<double> bath_lin;

    StatVector<double, SWE::n_variables> q_at_baryctr;
    HybMatrix<double, SWE::n_variables> q_at_vrtx;
    HybMatrix<double, SWE::n_variables> q_at_midpts;

    double bath_at_baryctr;
    DynRowVector<double> bath_at_vrtx;
    DynRowVector<double> bath_at_midpts;

    std::vector<bool> wet_neigh;
    AlignedVector<StatVector<double, SWE::n_variables>> q_at_baryctr_neigh;
    DynRowVector<double> bath_at_baryctr_neigh;

    HybMatrix<double, SWE::n_variables> delta;
    DynRowVector<double> bath_delta;

    AlignedVector<StatMatrix<double, SWE::n_dimensions, SWE::n_variables>> dq_r;
    HybMatrix<double, SWE::n_dimensions> dbath_r;
};
}

#endif
