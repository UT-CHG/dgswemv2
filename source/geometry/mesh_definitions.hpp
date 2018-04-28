#ifndef MESH_DEFINITIONS_HPP
#define MESH_DEFINITIONS_HPP

#include "mesh.hpp"
#include "raw_boundary.hpp"
#include "element.hpp"
#include "interface.hpp"
#include "boundary.hpp"

#include "../master/master_elements_2D.hpp"
#include "../shape/shapes_2D.hpp"
#include "../basis/bases_2D.hpp"
#include "../integration/integrations_1D.hpp"
#include "../integration/integrations_2D.hpp"

namespace Geometry {
template <typename Data>
using ElementTypeTuple = std::tuple<
    Element<2, Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>, Shape::StraightTriangle, Data>>;

template <typename Data, typename... ISs>
using InterfaceTypeTuple = std::tuple<Interface<1, Integration::GaussLegendre_1D, Data, ISs>...>;

template <typename Data, typename... BCs>
using BoundaryTypeTuple = std::tuple<Boundary<1, Integration::GaussLegendre_1D, Data, BCs>...>;

template <typename Data, typename Distributed>
using DistributedBoundaryTypeTuple = std::tuple<Boundary<1, Integration::GaussLegendre_1D, Data, Distributed>>;

template <typename Data, typename Distributed, typename IS, typename... BCs>
using MeshType = Mesh<ElementTypeTuple<Data>,
                      InterfaceTypeTuple<Data, IS>,
                      BoundaryTypeTuple<Data, BCs...>,
                      DistributedBoundaryTypeTuple<Data, Distributed>>;
};

#endif
