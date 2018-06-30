#ifndef MESH_SKELETON_HPP
#define MESH_SKELETON_HPP

#include "general_definitions.hpp"

#include "utilities//heterogeneous_containers.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename EdgeInterfaceTypeTuple, typename EdgeBoundaryTypeTuple>
class MeshSkeleton;

template <typename... EngeInternals, typename... EdgeBoundaries>
class MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>> {
  private:
    using EdgeInterfaceContainer = Utilities::HeterogeneousVector<EngeInternals...>;
    using EdgeBoundaryContainer  = Utilities::HeterogeneousVector<EdgeBoundaries...>;

    EdgeInterfaceContainer edge_interfaces;
    EdgeBoundaryContainer edge_boundaries;

  public:
    uint GetNumberEdgeInterfaces() { return this->edge_interfaces.size(); }
    uint GetNumberEdgeBoundaries() { return this->edge_boundaries.size(); }

    template <typename EdgeInterfaceType, typename... Args>
    void CreateEdgeInterface(Args&&... args);
    template <typename EdgeBoundaryType, typename... Args>
    void CreateEdgeBoundary(Args&&... args);

    template <typename F>
    void CallForEachEdgeInterface(const F& f);
    template <typename F>
    void CallForEachEdgeBoundary(const F& f);
};

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename EdgeInterfaceType, typename... Args>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CreateEdgeInterface(Args&&... args) {
    this->edge_interfaces.template emplace_back<EdgeInterfaceType>(std::forward<Args>(args)...);
}

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename EdgeBoundaryType, typename... Args>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CreateEdgeBoundary(Args&&... args) {
    this->edge_boundaries.template emplace_back<EdgeBoundaryType>(std::forward<Args>(args)...);
}

template <typename... EngeInternals, typename... EdgeBoundaries>
template <typename F>
void MeshSkeleton<std::tuple<EngeInternals...>, std::tuple<EdgeBoundaries...>>::CallForEachEdgeInterface(const F& f) {
    Utilities::for_each_in_tuple(this->edge_interfaces.data, [&f](auto& edge_interface_vector) {
        std::for_each(edge_interface_vector.begin(), edge_interface_vector.end(), f);
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