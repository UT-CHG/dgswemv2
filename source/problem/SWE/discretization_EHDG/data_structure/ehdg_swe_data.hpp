#ifndef EHDG_SWE_DATA_HPP
#define EHDG_SWE_DATA_HPP

#include "ehdg_swe_data_state.hpp"
#include "ehdg_swe_data_internal.hpp"
#include "ehdg_swe_data_boundary.hpp"
#include "ehdg_swe_data_source.hpp"
#include "ehdg_swe_data_wet_dry.hpp"
#include "ehdg_swe_data_slope_limit.hpp"

namespace SWE {
namespace EHDG {
struct Data {
    AlignedVector<State> state;
    Internal internal;
    AlignedVector<Boundary> boundary;

    Source source;
    WetDry wet_dry_state;
    SlopeLimit slope_limit_state;

    void initialize() {
        this->state.emplace_back(this->ndof);

        this->internal = Internal(this->ngp_internal);

        for (uint bound_id = 0; bound_id < this->nbound; ++bound_id) {
            this->boundary.push_back(Boundary(this->ngp_boundary[bound_id]));
        }

        this->source = Source(this->nnode);

        this->wet_dry_state = WetDry(this->nvrtx);

        this->slope_limit_state = SlopeLimit(this->nvrtx, this->nbound);
    }

    void resize(const uint nstate) {
        if (this->state.size() < nstate) {
            this->state.insert(this->state.end(), nstate - this->state.size(), State(ndof));
        } else if (this->state.size() > nstate) {
            this->state.erase(this->state.end() - (this->state.size() - nstate), this->state.end());
        }
    }

    uint get_nnode() { return this->nnode; }
    uint get_nvrtx() { return this->nvrtx; }
    uint get_nbound() { return this->nbound; }
    uint get_ndof() { return this->ndof; }
    uint get_ngp_internal() { return this->ngp_internal; }
    uint get_ngp_boundary(uint nbound) { return this->ngp_boundary[nbound]; }

    void set_nnode(const uint nnode) { this->nnode = nnode; }
    void set_nvrtx(const uint nvrtx) { this->nvrtx = nvrtx; }
    void set_nbound(const uint nbound) {
        this->nbound       = nbound;
        this->ngp_boundary = std::vector<uint>(this->nbound, 0);
    }
    void set_ndof(const uint ndof) { this->ndof = ndof; }
    void set_ngp_internal(const uint ngp) { this->ngp_internal = ngp; }
    void set_ngp_boundary(const uint bound_id, const uint ngp) { this->ngp_boundary[bound_id] = ngp; }

  private:
    uint nnode;
    uint nvrtx;
    uint nbound;
    uint ndof;
    uint ngp_internal;
    std::vector<uint> ngp_boundary;
};
}
}

#endif
