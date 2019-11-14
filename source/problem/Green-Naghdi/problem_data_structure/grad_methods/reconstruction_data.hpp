#ifndef RECONSTRUCTION_DATA_HPP
#define RECONSTRUCTION_DATA_HPP

namespace GN {
struct Reconstruction {
    Reconstruction() = default;
    Reconstruction(const uint nvrtx, const uint nbound, const std::vector<uint>& ngp_boundary)
        : local_nodeID(nvrtx),
          node_mult(nvrtx),
          dX_transpose(2, nbound),
          P(2, nbound),
          midpts_coord(nbound),
          baryctr_coord_neigh(nbound),
          a(nbound),
          a_elem(2 * nbound),
          /* *** */
          dze_at_midpts(GN::n_dimensions, nbound),
          dze_at_baryctr_neigh(nbound),
          dze_lin(GN::n_dimensions, nvrtx),
          /* *** */
          du_at_midpts(GN::n_du_terms, nbound),
          du_at_baryctr_neigh(nbound),
          du_lin(GN::n_du_terms, nvrtx),
          /* *** */
          ddu_at_midpts(GN::n_ddu_terms, nbound),
          ddu_at_baryctr_neigh(nbound),
          ddu_lin(GN::n_ddu_terms, nvrtx),
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
          dddbath_lin(GN::n_dddbath_terms, nvrtx) {
        // *** //
        set_constant(this->T, 1.0);
        this->T(0, 0) = -1.0;
        this->T(1, 1) = -1.0;
        this->T(2, 2) = -1.0;
    }

    std::vector<uint> local_nodeID;
    std::vector<uint> node_mult;

    HybMatrix<double, 2> dX_transpose;
    HybMatrix<double, 2> P;

    StatMatrix<double, SWE::n_variables, SWE::n_variables> T;

    double area;
    Point<2> baryctr_coord;
    AlignedVector<Point<2>> midpts_coord;
    AlignedVector<Point<2>> baryctr_coord_neigh;
    AlignedVector<StatVector<double, SWE::n_dimensions>> a;
    std::vector<uint> a_elem;

    StatVector<double, GN::n_dimensions> dze_at_baryctr;
    HybMatrix<double, GN::n_dimensions> dze_at_midpts;
    AlignedVector<StatVector<double, GN::n_dimensions>> dze_at_baryctr_neigh;
    HybMatrix<double, GN::n_dimensions> dze_lin;

    StatVector<double, GN::n_du_terms> du_at_baryctr;
    HybMatrix<double, GN::n_du_terms> du_at_midpts;
    AlignedVector<StatVector<double, GN::n_du_terms>> du_at_baryctr_neigh;
    HybMatrix<double, GN::n_du_terms> du_lin;

    StatVector<double, GN::n_ddu_terms> ddu_at_baryctr;
    HybMatrix<double, GN::n_ddu_terms> ddu_at_midpts;
    AlignedVector<StatVector<double, GN::n_ddu_terms>> ddu_at_baryctr_neigh;
    HybMatrix<double, GN::n_ddu_terms> ddu_lin;

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
};
}

#endif
