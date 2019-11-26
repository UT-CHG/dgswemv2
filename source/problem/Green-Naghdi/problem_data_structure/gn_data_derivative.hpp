#ifndef GN_DATA_DERIVATIVE_HPP
#define GN_DATA_DERIVATIVE_HPP

namespace GN {
struct Derivative {
    Derivative() = default;
    Derivative(const uint nvrtx, const uint nbound, const std::vector<uint>& ngp_boundary)
        : bath_hat_at_gp(nbound),
          ddbath_hat_at_gp(nbound),
          ze_hat_at_gp(nbound),
          u_hat_at_gp(nbound),
          du_hat_at_gp(nbound),
          /* *** */
          midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          /* *** */
          a(nbound),
          a_elem(2 * nbound),
          /* *** */
          dX_transpose(2, nbound),
          P(2, nbound),
          /* *** */
          local_nodeID(nvrtx),
          node_mult(nvrtx),
          /* *** */
          bath_at_midpts(nbound),
          bath_at_baryctr_neigh(nbound),
          bath_lin(nvrtx),
          /* *** */
          dbath_at_midpts(GN::n_dimensions, nbound),
          dbath_at_baryctr_neigh(nbound),
          dbath_lin(GN::n_dimensions, nvrtx),
          /* *** */
          ddbath_at_midpts(GN::n_ddbath_terms, nbound),
          ddbath_at_baryctr_neigh(nbound),
          ddbath_lin(GN::n_ddbath_terms, nvrtx),
          /* *** */
          dddbath_at_midpts(GN::n_dddbath_terms, nbound),
          dddbath_at_baryctr_neigh(nbound),
          dddbath_lin(GN::n_dddbath_terms, nvrtx),
          /* *** */
          ze_at_midpts(nbound),
          ze_at_baryctr_neigh(nbound),
          ze_lin(nvrtx),
          /* *** */
          dze_at_midpts(GN::n_dimensions, nbound),
          dze_at_baryctr_neigh(nbound),
          dze_lin(GN::n_dimensions, nvrtx),
          /* *** */
          u_at_midpts(GN::n_dimensions, nbound),
          u_at_baryctr_neigh(nbound),
          u_lin(GN::n_dimensions, nvrtx),
          /* *** */
          du_at_midpts(GN::n_du_terms, nbound),
          du_at_baryctr_neigh(nbound),
          du_lin(GN::n_du_terms, nvrtx),
          /* *** */
          ddu_at_midpts(GN::n_ddu_terms, nbound),
          ddu_at_baryctr_neigh(nbound),
          ddu_lin(GN::n_ddu_terms, nvrtx) {
        /* *** */
        for (uint bound = 0; bound < nbound; ++bound) {
            bath_hat_at_gp[bound]   = DynRowVector<double>(ngp_boundary[bound]);
            ddbath_hat_at_gp[bound] = HybMatrix<double, GN::n_ddbath_terms>(GN::n_ddbath_terms, ngp_boundary[bound]);

            bath_hat_at_gp[bound] = DynRowVector<double>(ngp_boundary[bound]);
            u_hat_at_gp[bound]    = HybMatrix<double, GN::n_dimensions>(GN::n_dimensions, ngp_boundary[bound]);
            du_hat_at_gp[bound]   = HybMatrix<double, GN::n_du_terms>(GN::n_du_terms, ngp_boundary[bound]);
        }

        set_constant(this->T, 1.0);
        this->T(0, 0) = -1.0;
        this->T(1, 1) = -1.0;
        this->T(2, 2) = -1.0;
    }

    /* L2 projection data for LDG */
    std::vector<DynRowVector<double>> bath_hat_at_gp;
    AlignedVector<HybMatrix<double, GN::n_ddbath_terms>> ddbath_hat_at_gp;

    std::vector<DynRowVector<double>> ze_hat_at_gp;
    AlignedVector<HybMatrix<double, GN::n_dimensions>> u_hat_at_gp;
    AlignedVector<HybMatrix<double, GN::n_du_terms>> du_hat_at_gp;

    /* Geometry */
    double area;
    Point<2> baryctr_coord;
    AlignedVector<Point<2>> midpts_coord;
    AlignedVector<Point<2>> baryctr_coord_neigh;

    /* Linear interpolation data */
    AlignedVector<StatVector<double, SWE::n_dimensions>> a;
    std::vector<uint> a_elem;

    /* Least squares data */
    HybMatrix<double, 2> dX_transpose;
    HybMatrix<double, 2> P;

    /* Averaging data */
    std::vector<uint> local_nodeID;
    std::vector<uint> node_mult;

    /* Transform from midpoints to vertexes */
    StatMatrix<double, SWE::n_variables, SWE::n_variables> T;

    double bath_at_baryctr;
    DynRowVector<double> bath_at_midpts;
    DynVector<double> bath_at_baryctr_neigh;
    DynRowVector<double> bath_lin;

    StatVector<double, GN::n_dimensions> dbath_at_baryctr;
    HybMatrix<double, GN::n_dimensions> dbath_at_midpts;
    AlignedVector<StatVector<double, GN::n_dimensions>> dbath_at_baryctr_neigh;
    HybMatrix<double, GN::n_dimensions> dbath_lin;

    StatVector<double, GN::n_ddbath_terms> ddbath_at_baryctr;
    HybMatrix<double, GN::n_ddbath_terms> ddbath_at_midpts;
    AlignedVector<StatVector<double, GN::n_ddbath_terms>> ddbath_at_baryctr_neigh;
    HybMatrix<double, GN::n_ddbath_terms> ddbath_lin;

    StatVector<double, GN::n_dddbath_terms> dddbath_at_baryctr;
    HybMatrix<double, GN::n_dddbath_terms> dddbath_at_midpts;
    AlignedVector<StatVector<double, GN::n_dddbath_terms>> dddbath_at_baryctr_neigh;
    HybMatrix<double, GN::n_dddbath_terms> dddbath_lin;

    double ze_at_baryctr;
    DynRowVector<double> ze_at_midpts;
    DynVector<double> ze_at_baryctr_neigh;
    DynRowVector<double> ze_lin;

    StatVector<double, GN::n_dimensions> dze_at_baryctr;
    HybMatrix<double, GN::n_dimensions> dze_at_midpts;
    AlignedVector<StatVector<double, GN::n_dimensions>> dze_at_baryctr_neigh;
    HybMatrix<double, GN::n_dimensions> dze_lin;

    StatVector<double, GN::n_dimensions> u_at_baryctr;
    HybMatrix<double, GN::n_dimensions> u_at_midpts;
    AlignedVector<StatVector<double, GN::n_dimensions>> u_at_baryctr_neigh;
    HybMatrix<double, GN::n_dimensions> u_lin;

    StatVector<double, GN::n_du_terms> du_at_baryctr;
    HybMatrix<double, GN::n_du_terms> du_at_midpts;
    AlignedVector<StatVector<double, GN::n_du_terms>> du_at_baryctr_neigh;
    HybMatrix<double, GN::n_du_terms> du_lin;

    StatVector<double, GN::n_ddu_terms> ddu_at_baryctr;
    HybMatrix<double, GN::n_ddu_terms> ddu_at_midpts;
    AlignedVector<StatVector<double, GN::n_ddu_terms>> ddu_at_baryctr_neigh;
    HybMatrix<double, GN::n_ddu_terms> ddu_lin;
};
}

#endif
