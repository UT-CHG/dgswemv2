#ifndef MESH_DEFINITIONS_H
#define MESH_DEFINITIONS_H

#include "general_definitions.hpp"

#include "class_element.hpp"
#include "class_boundary.hpp"

namespace Geometry {
	template<typename Data>
	using ElementTypeTuple = std::tuple<
		Element<2, Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>, Shape::StraightTriangle, Data>
	>;

	template<typename Data>
	using InterfaceTypeTuple = std::tuple<
		Interface<1, Integration::GaussLegendre_1D, Data>
	>;

	template<typename Data>
	using BoundaryTypeTuple = std::tuple<
		Boundary<1, Integration::GaussLegendre_1D, Data>
	>;

	template<typename Data>
	using MeshType = Mesh<ElementTypeTuple<Data>, InterfaceTypeTuple<Data>, BoundaryTypeTuple<Data>>;
};

#endif