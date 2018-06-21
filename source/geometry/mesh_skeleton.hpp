#ifndef MESH_SKELETON_HPP
#define MESH_SKELETON_HPP

#include "general_definitions.hpp"

#include "utilities//heterogeneous_containers.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename EdgeInternalTypeTuple, typename EdgeBoundaryTypeTuple>
class MeshSkeleton;

template <typename... EngeInternals, typename... EdgeBoundaries>
class MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>> {
  private:
    using EdgeInternalContainer = Utilities::HeterogeneousVector<EngeInternals...>;
    using EdgeBoundaryContainer = Utilities::HeterogeneousVector<EdgeBoundaries...>;

    EdgeInternalContainer edge_internals;
    EdgeBoundaryContainer edge_boundaries;

  public:
    uint GetNumberEdgeInternals() { return this->edge_internals.size(); }
    uint GetNumberEdgeBoundaries() { return this->edge_boundaries.size(); }

    template <typename EdgeInternalType, typename... Args>
    void CreateEdgeInternal(Args&&... args);
    template <typename EdgeBoundaryType, typename... Args>
    void CreateEdgeBoundary(Args&&... args);

    template <typename F>
    void CallForEachEdgeInternal(const F& f);
    template <typename F>
    void CallForEachEdgeBoundary(const F& f);
};

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename EdgeInternalType, typename... Args>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CreateEdgeInternal(Args&&... args) {
    this->edge_internals.template emplace_back<EdgeInternalType>(std::forward<Args>(args)...);
}

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename EdgeBoundaryType, typename... Args>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CreateEdgeBoundary(Args&&... args) {
    this->edge_boundaries.template emplace_back<EdgeBoundaryType>(std::forward<Args>(args)...);
}

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename F>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CallForEachEdgeInternal(const F& f) {
    Utilities::for_each_in_tuple(this->edge_internals.data, [&f](auto& edge_internal_vector) {
        std::for_each(edge_internal_vector.begin(), edge_internal_vector.end(), f);
    });
}

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename F>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CallForEachEdgeBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->edge_boundaries.data, [&f](auto& edge_boundary_vector) {
        std::for_each(edge_boundary_vector.begin(), edge_boundary_vector.end(), f);
    });
}
}

#endif