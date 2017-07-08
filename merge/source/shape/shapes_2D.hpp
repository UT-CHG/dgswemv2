#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "../general_definitions.hpp"

namespace Shape {
	class StraightTriangle : public Shape<2> {
	public:
		StraightTriangle(std::vector<Point<2>> nodal_coordinates) : Shape<2>(nodal_coordinates) {}

		std::vector<double> GetJdet(const std::vector<Point<2>>&);
		Array3D<double> GetJinv(const std::vector<Point<2>>&);
		std::vector<double> GetSurfaceJ(uint, const std::vector<Point<2>>&);
		Array2D<double> GetSurfaceNormal(uint, const std::vector<Point<2>>&);
		void GetVTK(std::vector<Point<3>>&, Array2D<uint>&);
	};
}

#endif