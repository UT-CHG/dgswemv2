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
        set_constant(this->T, 1.0);
        this->T(0, 0) = -1.0;
        this->T(1, 1) = -1.0;
        this->T(2, 2) = -1.0;
    }

    std::vector<double> lengths;
    double radius;
    bool troubled    = false;
    double perimeter = 0.0;
    double I         = 0.0;
    std::vector<DynRowVector<double>> hdif_at_gp;

    StatMatrix<double, SWE::n_variables, SWE::n_variables> T;

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

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & T
            & baryctr_coord
            & midpts_coord
            & baryctr_coord_neigh
            & median
            & alpha
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

#endif
