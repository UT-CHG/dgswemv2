#ifndef SWE_DATA_HPP
#define SWE_DATA_HPP

#include "../../general_definitions.hpp"

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
};

struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : ze_flux_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          qx_flux_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          qy_flux_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          ze_source_term_at_gp(ngp),
          qx_source_term_at_gp(ngp),
          qy_source_term_at_gp(ngp),
          ze_at_gp(ngp),
          qx_at_gp(ngp),
          qy_at_gp(ngp),
          bath_at_gp(ngp),
          bath_deriv_wrt_x_at_gp(ngp),
          bath_deriv_wrt_y_at_gp(ngp),
          h_at_gp(ngp) {}

    std::array<std::vector<double>, 2> ze_flux_at_gp;
    std::array<std::vector<double>, 2> qx_flux_at_gp;
    std::array<std::vector<double>, 2> qy_flux_at_gp;

    std::vector<double> ze_source_term_at_gp;
    std::vector<double> qx_source_term_at_gp;
    std::vector<double> qy_source_term_at_gp;

    std::vector<double> ze_at_gp;
    std::vector<double> qx_at_gp;
    std::vector<double> qy_at_gp;

    std::vector<double> bath_at_gp;
    std::vector<double> bath_deriv_wrt_x_at_gp;
    std::vector<double> bath_deriv_wrt_y_at_gp;

    std::vector<double> h_at_gp;
};

struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : ze_at_gp(ngp),
          qx_at_gp(ngp),
          qy_at_gp(ngp),
          bath_at_gp(ngp),
          ze_numerical_flux_at_gp(ngp),
          qx_numerical_flux_at_gp(ngp),
          qy_numerical_flux_at_gp(ngp) {}

    std::vector<double> ze_at_gp;
    std::vector<double> qx_at_gp;
    std::vector<double> qy_at_gp;

    std::vector<double> bath_at_gp;

    std::vector<double> ze_numerical_flux_at_gp;
    std::vector<double> qx_numerical_flux_at_gp;
    std::vector<double> qy_numerical_flux_at_gp;
};

struct WetDry {
    WetDry() = default;
    WetDry(const uint nvrtx)
        : ze_at_vrtx(nvrtx),
          qx_at_vrtx(nvrtx),
          qy_at_vrtx(nvrtx),
          bath_at_vrtx(nvrtx),
          h_at_vrtx(nvrtx),
          h_at_vrtx_temp(nvrtx) {}

    bool wet;

    double bath_min;

    std::vector<double> ze_at_vrtx;
    std::vector<double> qx_at_vrtx;
    std::vector<double> qy_at_vrtx;

    std::vector<double> bath_at_vrtx;

    std::vector<double> h_at_vrtx;
    std::vector<double> h_at_vrtx_temp;
};

struct Data {
    WetDry wet_dry_state;

    std::vector<State> state;
    Internal internal;
    std::vector<Boundary> boundary;

    void initialize() {
        this->wet_dry_state = WetDry(this->nvrtx);

        this->state = std::vector<State>{State(this->ndof)};

        this->internal = Internal(this->ngp_internal);

        for (uint bound_id = 0; bound_id < this->nbound; bound_id++) {
            this->boundary.push_back(Boundary(this->ngp_boundary[bound_id]));
        }
    }

    void resize(const uint nstate) {
        if ((this->state.size() - 1) < nstate) {
            this->state.insert(this->state.end(), nstate - (this->state.size() - 1), State(ndof));
        } else if ((this->state.size() - 1) > nstate) {
            this->state.erase(this->state.end() - ((this->state.size() - 1) - nstate), this->state.end());
        }
    }

    uint get_nvrtx() { return this->nvrtx; }
    uint get_ndof() { return this->ndof; }
    uint get_ngp_internal() { return this->ngp_internal; }
    uint get_nbound() { return this->nbound; }
    uint get_ngp_boundary(uint nbound) { return this->ngp_boundary[nbound]; }

    void set_nvrtx(const uint nvrtx) { this->nvrtx = nvrtx; }
    void set_ndof(const uint ndof) { this->ndof = ndof; }
    void set_ngp_internal(const uint ngp) { this->ngp_internal = ngp; }
    void set_nbound(const uint nbound) {
        this->nbound = nbound;
        this->ngp_boundary = std::vector<uint>(this->nbound, 0);
    }
    void set_ngp_boundary(const uint bound_id, const uint ngp) { this->ngp_boundary[bound_id] = ngp; }

  private:
    uint nvrtx;
    uint ndof;
    uint ngp_internal;
    uint nbound;
    std::vector<uint> ngp_boundary;
};
}

#endif
