#ifndef MESH_SKELETON_HPP
#define MESH_SKELETON_HPP

#include "general_definitions.hpp"

#include "utilities//heterogeneous_containers.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename EdgeInterfaceTypeTuple, typename EdgeBoundaryTypeTuple, typename EdgeDistributedTypeTuple>
class MeshSkeleton;

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
class MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>> {
  private:
    using EdgeInterfaceContainer   = Utilities::HeterogeneousVector<EdgeInterfaces...>;
    using EdgeBoundaryContainer    = Utilities::HeterogeneousVector<EdgeBoundaries...>;
    using EdgeDistributedContainer = Utilities::HeterogeneousVector<EdgeDistributeds...>;

    EdgeInterfaceContainer edge_interfaces;
    EdgeBoundaryContainer edge_boundaries;
    EdgeDistributedContainer edge_distributeds;

  public:
    uint GetNumberEdgeInterfaces() { return this->edge_interfaces.size(); }
    uint GetNumberEdgeBoundaries() { return this->edge_boundaries.size(); }
    uint GetNumberEdgeDistributeds() { return this->edge_distributeds.size(); }

    template <typename EdgeInterfaceType, typename... Args>
    void CreateEdgeInterface(Args&&... args);
    template <typename EdgeBoundaryType, typename... Args>
    void CreateEdgeBoundary(Args&&... args);
    template <typename EdgeDistributedType, typename... Args>
    void CreateEdgeDistributed(Args&&... args);

    template <typename F>
    void CallForEachEdgeInterface(const F& f);
    template <typename F>
    void CallForEachEdgeBoundary(const F& f);
    template <typename F>
    void CallForEachEdgeDistributed(const F& f);
};

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
template <typename EdgeInterfaceType, typename... Args>
void MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>>::
    CreateEdgeInterface(Args&&... args) {
    this->edge_interfaces.template emplace_back<EdgeInterfaceType>(std::forward<Args>(args)...);
}

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
template <typename EdgeBoundaryType, typename... Args>
void MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>>::
    CreateEdgeBoundary(Args&&... args) {
    this->edge_boundaries.template emplace_back<EdgeBoundaryType>(std::forward<Args>(args)...);
}

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
template <typename EdgeDistributedType, typename... Args>
void MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>>::
    CreateEdgeDistributed(Args&&... args) {
    this->edge_distributeds.template emplace_back<EdgeDistributedType>(std::forward<Args>(args)...);
}

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
template <typename F>
void MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>>::
    CallForEachEdgeInterface(const F& f) {
    Utilities::for_each_in_tuple(this->edge_interfaces.data, [&f](auto& edge_interface_vector) {
        std::for_each(edge_interface_vector.begin(), edge_interface_vector.end(), f);
    });
}

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
template <typename F>
void MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>>::
    CallForEachEdgeBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->edge_boundaries.data, [&f](auto& edge_boundary_vector) {
        std::for_each(edge_boundary_vector.begin(), edge_boundary_vector.end(), f);
    });
}

template <typename... EdgeInterfaces, typename... EdgeBoundaries, typename... EdgeDistributeds>
template <typename F>
void MeshSkeleton<std::tuple<EdgeInterfaces...>, std::tuple<EdgeBoundaries...>, std::tuple<EdgeDistributeds...>>::
    CallForEachEdgeDistributed(const F& f) {
    Utilities::for_each_in_tuple(this->edge_distributeds.data, [&f](auto& edge_distributed_vector) {
        std::for_each(edge_distributed_vector.begin(), edge_distributed_vector.end(), f);
    });
}
}

#endif