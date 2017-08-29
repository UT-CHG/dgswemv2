#ifndef SWE_DATA_HPP
#define SWE_DATA_HPP

#include <array>
#include <vector>

namespace SWE {
struct State {
    State() = default;
    State(uint ndof) : ze(ndof), qx(ndof), qy(ndof), bath(ndof), rhs_ze(ndof), rhs_qx(ndof), rhs_qy(ndof) {}

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
    Internal(uint ngp)
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
          water_column_hgt_at_gp(ngp) {}

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

    std::vector<double> water_column_hgt_at_gp;
};

struct Boundary {
    Boundary() = default;
    Boundary(uint ngp)
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

struct Data {
    void initialize() {
        this->state = std::vector<State>{State(this->ndof)};

        this->internal = Internal(this->ngp_internal);

        for (uint bound = 0; bound < this->nbound; bound++) {
            this->boundary.push_back(Boundary(this->ngp_boundary[bound]));
        }
    }

    void resize(uint nstate) {
        if ((this->state.size() - 1) < nstate) {
            this->state.insert(this->state.end(), nstate - (this->state.size() - 1), State(ndof));
        } else if ((this->state.size() - 1) > nstate) {
            this->state.erase(this->state.end() - ((this->state.size() - 1) - nstate), this->state.end());
        }
    }

    std::vector<State> state;
    Internal internal;
    std::vector<Boundary> boundary;

    uint get_ndof() { return this->ndof; }
    uint get_ngp_internal() { return this->ngp_internal; }
    uint get_nbound() { return this->nbound; }
    uint get_ngp_boundary(uint nbound) { return this->ngp_boundary[nbound]; }

    void set_ndof(uint ndof) { this->ndof = ndof; }
    void set_ngp_internal(uint ngp) { this->ngp_internal = ngp; }
    void set_nbound(uint nbound) {
        this->nbound = nbound;
        this->ngp_boundary = std::vector<uint>(this->nbound, 0);
    }
    void set_ngp_boundary(uint nbound, uint ngp) { this->ngp_boundary[nbound] = ngp; }

  private:
    uint ndof;
    uint ngp_internal;
    uint nbound;
    std::vector<uint> ngp_boundary;
};
}

#endif
