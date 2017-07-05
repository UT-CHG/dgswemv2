#ifndef CLASS_MASTER_ELEMENT_HPP
#define CLASS_MASTER_ELEMENT_HPP

#include "../general_definitions.hpp"

#include "../basis/bases_2D.hpp"

#include "../integration/integrations_1D.hpp"
#include "../integration/integrations_2D.hpp"

namespace Master {
	template<class basis_type, class integration_type>
	class Triangle : public Master<2> {
	public:
		basis_type basis;
		integration_type integration;

	public:
		Triangle(uint);

		std::vector<Point<2>> BoundaryToMasterCoordinates(uint, const std::vector<Point<1>>&);

	private:
		std::vector<Point<2>> VTKPostCell();
		std::vector<Point<2>> VTKPostPoint();
	};
}

#include "elements_2D/master_triangle.tpp"

#endif