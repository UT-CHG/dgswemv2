#ifndef RKDG_SWE_DATA_HPP
#define RKDG_SWE_DATA_HPP

#include "general_definitions.hpp"

#include "rkdg_swe_data_state.hpp"
#include "rkdg_swe_data_internal.hpp"
#include "rkdg_swe_data_boundary.hpp"
#include "rkdg_swe_data_source.hpp"
#include "rkdg_swe_data_wet_dry.hpp"
#include "rkdg_swe_data_slope_limit.hpp"

namespace SWE {
namespace RKDG {
struct Data {
    std::vector<State> state;
    Internal internal;
    std::vector<Boundary> boundary;

    Source source;
    WetDry wet_dry_state;
    SlopeLimit slope_limit_state;

    void initialize() {
        this->state = std::vector<State>{State(this->ndof)};

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

    uint get_nnode() const { return this->nnode; }
    uint get_nvrtx() const { return this->nvrtx; }
    uint get_nbound() const { return this->nbound; }
    uint get_ndof() const { return this->ndof; }
    uint get_ngp_internal() const { return this->ngp_internal; }
    uint get_ngp_boundary(uint nbound) const { return this->ngp_boundary[nbound]; }

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

  public:
#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & state
            & internal
            & boundary
            & source
            & wet_dry_state
            & slope_limit_state
            & nnode
            & nvrtx
            & nbound
            & ndof
            & ngp_internal
            & ngp_boundary;
        // clang-format on
    }
#endif
};
}
}

#endif
