#ifndef MESH_SKELETON_DEFINITIONS_HPP
#define MESH_SKELETON_DEFINITIONS_HPP

#include "general_definitions.hpp"

#include "mesh_skeleton.hpp"
#include "edge_interface.hpp"
#include "edge_boundary.hpp"

#include "basis/bases_1D.hpp"
#include "basis/bases_2D.hpp"
#include "integration/integrations_1D.hpp"
#include "integration/integrations_2D.hpp"

namespace Geometry {
template <typename EdgeData, typename InterfaceType>
struct EdgeInterfaceTypeTuple;

template <typename EdgeData, typename... InterfaceTypes>
struct EdgeInterfaceTypeTuple<EdgeData, std::tuple<InterfaceTypes...>> {
    using Type = std::tuple<EdgeInterface<1, Basis::Legendre_1D, EdgeData, InterfaceTypes>...>;
};

template <typename EdgeData, typename BoundaryType>
struct EdgeBoundaryTypeTuple;

template <typename EdgeData, typename... BoundaryTypes>
struct EdgeBoundaryTypeTuple<EdgeData, std::tuple<BoundaryTypes...>> {
    using Type = std::tuple<EdgeBoundary<1, Basis::Legendre_1D, EdgeData, BoundaryTypes>...>;
};

template <typename EdgeData, typename DistributedBoundaryType>
struct EdgeDistributedTypeTuple;

template <typename EdgeData, typename... DistributedBoundaryTypes>
struct EdgeDistributedTypeTuple<EdgeData, std::tuple<DistributedBoundaryTypes...>> {
    using Type = std::tuple<EdgeBoundary<1, Basis::Legendre_1D, EdgeData, DistributedBoundaryTypes>...>;
};

template <typename EdgeData, typename InterfaceType, typename BoundaryType, typename DistributedBoundaryType>
struct MeshSkeletonType;

template <typename EdgeData,
          typename... InterfaceTypes,
          typename... BoundaryTypes,
          typename... DistributedBoundaryTypes>
struct MeshSkeletonType<EdgeData,
                        std::tuple<InterfaceTypes...>,
                        std::tuple<BoundaryTypes...>,
                        std::tuple<DistributedBoundaryTypes...>> {
    using Type =
        MeshSkeleton<typename EdgeInterfaceTypeTuple<EdgeData, std::tuple<InterfaceTypes...>>::Type,
                     typename EdgeBoundaryTypeTuple<EdgeData, std::tuple<BoundaryTypes...>>::Type,
                     typename EdgeDistributedTypeTuple<EdgeData, std::tuple<DistributedBoundaryTypes...>>::Type>;
};
}

#endif
