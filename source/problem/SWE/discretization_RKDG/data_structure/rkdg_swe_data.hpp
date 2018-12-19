#ifndef RKDG_SWE_DATA_HPP
#define RKDG_SWE_DATA_HPP

#include "rkdg_swe_data_state.hpp"
#include "rkdg_swe_data_internal.hpp"
#include "rkdg_swe_data_boundary.hpp"
#include "rkdg_swe_data_source.hpp"
#include "rkdg_swe_data_wet_dry.hpp"
#include "rkdg_swe_data_slope_limit.hpp"

#include "utilities/at_each.hpp"

namespace SWE {
namespace RKDG {
struct Accessor {
    Accessor(AlignedVector<StateAccessor>&& state_,
             InternalAccessor&& internal_,
             AlignedVector<BoundaryAccessor>&& boundary_)
        : state(state_), internal(internal_), boundary(boundary_)
    {}

    AlignedVector<StateAccessor> state;
    InternalAccessor internal;
    AlignedVector<BoundaryAccessor> boundary;

    Source source;
    WetDry wet_dry_state;
    SlopeLimit slope_limit_state;

    void initialize() {
        this->source = Source(this->nnode);
        this->wet_dry_state = WetDry(this->nvrtx);

        this->slope_limit_state = SlopeLimit(this->nvrtx, this->nbound);
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

struct SoAContainer {
    using AccessorType = Accessor;

    SoAContainer()=default;
    SoAContainer(uint ndofs,
                 uint ngp_internal,
                 uint ngp_edge,
                 uint nstages,
                 uint nelements,
                 uint nsides);

    Accessor at(const size_t index);

    AlignedVector<StateData> state;
    InternalData internal;
    AlignedVector<BoundaryData> boundary;
};

SoAContainer::SoAContainer(uint ndofs,
                           uint ngp_internal,
                           uint ngp_edge,
                           uint nstages,
                           uint nelements,
                           uint nsides)
    : internal(nelements, ngp_internal) {
    this->state.reserve(nstages+1);
    for ( uint i = 0; i < nstages+1; ++i ) {
        this->state.emplace_back(nelements, ndofs);
    }

    for ( uint side = 0; side < nsides; ++side ) {
        this->boundary.emplace_back(nelements, ngp_edge);
    }
}

Accessor SoAContainer::at(const size_t index) {
    return Accessor(std::move(Utilities::at_each(this->state, index)),
                    std::move(this->internal.at(index)),
                    std::move(Utilities::at_each(this->boundary,index))
                   );
}

}
}

#endif
