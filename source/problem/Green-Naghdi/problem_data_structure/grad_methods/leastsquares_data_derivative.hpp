#ifndef LEASTSQUARES_DATA_DERIVATIVE_HPP
#define LEASTSQUARES_DATA_DERIVATIVE_HPP

#include "reconstruction_data.hpp"

namespace GN {
struct Derivative : GN::Reconstruction {
    Derivative() = default;
    Derivative(const uint nvrtx, const uint nbound, const std::vector<uint>& ngp_boundary)
        : GN::Reconstruction(nvrtx, nbound, ngp_boundary),
          bath_at_midpts(nbound),
          bath_at_baryctr_neigh(nbound),
          bath_lin(nvrtx),
          ze_at_midpts(nbound),
          ze_at_baryctr_neigh(nbound),
          ze_lin(nvrtx),
          u_at_midpts(GN::n_dimensions, nbound),
          u_at_baryctr_neigh(nbound),
          u_lin(GN::n_dimensions, nvrtx) {}

    double bath_at_baryctr;
    DynRowVector<double> bath_at_midpts;
    DynVector<double> bath_at_baryctr_neigh;
    DynRowVector<double> bath_lin;

    double ze_at_baryctr;
    DynRowVector<double> ze_at_midpts;
    DynVector<double> ze_at_baryctr_neigh;
    DynRowVector<double> ze_lin;

    StatVector<double, GN::n_dimensions> u_at_baryctr;
    HybMatrix<double, GN::n_dimensions> u_at_midpts;
    AlignedVector<StatVector<double, GN::n_dimensions>> u_at_baryctr_neigh;
    HybMatrix<double, GN::n_dimensions> u_lin;
};
}

#endif
