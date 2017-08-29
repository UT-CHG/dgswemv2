#ifndef MESH_DEFINITIONS_H
#define MESH_DEFINITIONS_H

#include "mesh.hpp"
#include "element.hpp"
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

template <typename Data>
using InterfaceTypeTuple = std::tuple<Interface<1, Integration::GaussLegendre_1D, Data>>;

template <typename Data, typename... BCs>
using BoundaryTypeTuple = std::tuple<Boundary<1, Integration::GaussLegendre_1D, Data, BCs>...>;

template <typename Data, typename Distributed>
using DistributedInterface = std::tuple<Boundary<1, Integration::GaussLegendre_1D, Data, Distributed>>;

template <typename Data, typename Distributed, typename... BCs>
using MeshType = Mesh<ElementTypeTuple<Data>,
                      InterfaceTypeTuple<Data>,
                      BoundaryTypeTuple<Data, BCs...>,
                      DistributedInterface<Data, Distributed>>;
};

#endif
