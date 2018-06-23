#ifndef MESH_SKELETON_DEFINITIONS_HPP
#define MESH_SKELETON_DEFINITIONS_HPP

#include "mesh_skeleton.hpp"
#include "edge_internal.hpp"
#include "edge_boundary.hpp"

#include "basis/bases_1D.hpp"
#include "basis/bases_2D.hpp"
#include "integration/integrations_1D.hpp"
#include "integration/integrations_2D.hpp"

namespace Geometry {
template <typename Data, typename EdgeData>
using EdgeInternalTypeTuple = std::tuple<EdgeInternal<1, Basis::Legendre_1D, Data, EdgeData>>;

template <typename Data, typename EdgeData>
using EdgeBoundaryTypeTuple = std::tuple<EdgeBoundary<1, Basis::Legendre_1D, Data, EdgeData>>;

template <typename Data, typename EdgeData>
using MeshSkeletonType = MeshSkeleton<EdgeInternalTypeTuple<Data, EdgeData>, EdgeBoundaryTypeTuple<Data, EdgeData>>;
}

#endif
