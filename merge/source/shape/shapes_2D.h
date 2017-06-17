#ifndef SHAPES_2D_H
#define SHAPES_2D_H

#include "../general_definitions.h"

namespace Shape{
class StraightTriangle : Shape<2> {
public:
	std::vector<double> get_J_det(const std::vector<Point<2>>&);
	Array3D<double> get_J_inv(const std::vector<Point<2>>&);
	Array2D<double> get_surface_J(const std::vector<Point<2>>&);
	Array3D<double> get_surface_normal(const std::vector<Point<2>>&);
};
}

#endif