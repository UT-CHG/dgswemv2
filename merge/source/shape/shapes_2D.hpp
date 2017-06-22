#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "../general_definitions.hpp"

namespace Shape {
	class StraightTriangle : Shape<2> {
	public:
		std::vector<double> get_J_det(const std::vector<Point<2>>&);
		Array3D<double> get_J_inv(const std::vector<Point<2>>&);
		std::vector<double> get_surface_J(uint, const std::vector<Point<2>>&);
		Array2D<double> get_surface_normal(uint, const std::vector<Point<2>>&);
		void get_VTK(std::vector<Point<3>>&, Array2D<uint>&, const std::vector<Point<2>>&);
	};
}

#endif