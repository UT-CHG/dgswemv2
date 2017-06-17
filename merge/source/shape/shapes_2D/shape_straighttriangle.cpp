#include "../shapes_2D.h"

namespace Shape{
std::vector<double> StraightTriangle::get_J_det(const std::vector<Point<2>>& nodal_coordinates){
		std::vector<double> J_det;
        
        Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((nodal_coordinates[1][X] - nodal_coordinates[0][X]) / 2.0);
		J[0].push_back((nodal_coordinates[2][X] - nodal_coordinates[0][X]) / 2.0);
		J[1].push_back((nodal_coordinates[1][Y] - nodal_coordinates[0][Y]) / 2.0);
		J[1].push_back((nodal_coordinates[2][Y] - nodal_coordinates[0][Y]) / 2.0);

		J_det.push_back(J[0][0] * J[1][1] - J[0][1] * J[1][0]);

        return J_det;
}

Array3D<double> StraightTriangle::get_J_inv(const std::vector<Point<2>>& nodal_coordinates){
		Array3D<double> J_inv(2);
        J_inv[0].reserve(2);
        J_inv[1].reserve(2);

        Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((nodal_coordinates[1][X] - nodal_coordinates[0][X]) / 2.0);
		J[0].push_back((nodal_coordinates[2][X] - nodal_coordinates[0][X]) / 2.0);
		J[1].push_back((nodal_coordinates[1][Y] - nodal_coordinates[0][Y]) / 2.0);
		J[1].push_back((nodal_coordinates[2][Y] - nodal_coordinates[0][Y]) / 2.0);

		double det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0];

		J_inv[0][0].push_back(J[1][1] / det_J);
		J_inv[0][1].push_back(-J[0][1] / det_J);
		J_inv[1][0].push_back(-J[1][0] / det_J);
		J_inv[1][1].push_back(J[0][0] / det_J);

        return J_inv;
}

Array2D<double> StraightTriangle::get_surface_J(const std::vector<Point<2>>& nodal_coordinates){
		Array2D<double> surface_J(3); 
        
        Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((nodal_coordinates[1][X] - nodal_coordinates[0][X]) / 2.0);
		J[0].push_back((nodal_coordinates[2][X] - nodal_coordinates[0][X]) / 2.0);
		J[1].push_back((nodal_coordinates[1][Y] - nodal_coordinates[0][Y]) / 2.0);
		J[1].push_back((nodal_coordinates[2][Y] - nodal_coordinates[0][Y]) / 2.0);

		surface_J[0].push_back(sqrt(pow((J[0][0] - J[0][1]), 2.0) + pow((J[1][0] - J[1][1]), 2.0)));
		surface_J[1].push_back(sqrt(pow(J[0][1], 2.0) + pow(J[1][1], 2.0)));
		surface_J[2].push_back(sqrt(pow(J[0][0], 2.0) + pow(J[1][0], 2.0)));

        return surface_J;
}

Array3D<double> StraightTriangle::get_surface_normal(const std::vector<Point<2>>& nodal_coordinates){
		Array3D<double> surface_normal(3);
        surface_normal[0].reserve(2);
        surface_normal[1].reserve(2);
        surface_normal[2].reserve(2);
        
        Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((nodal_coordinates[1][X] - nodal_coordinates[0][X]) / 2.0);
		J[0].push_back((nodal_coordinates[2][X] - nodal_coordinates[0][X]) / 2.0);
		J[1].push_back((nodal_coordinates[1][Y] - nodal_coordinates[0][Y]) / 2.0);
		J[1].push_back((nodal_coordinates[2][Y] - nodal_coordinates[0][Y]) / 2.0);

   		Array2D<double> surface_J = this->get_surface_J(nodal_coordinates); 

		double det_J = this->get_J_det(nodal_coordinates)[0];

		double cw = det_J / abs(det_J); //CW or CCW

		surface_normal[0][X].push_back((J[1][1] - J[1][0]) / surface_J[0][0] * cw);
		surface_normal[1][X].push_back(-J[1][1] / surface_J[1][0] * cw);
		surface_normal[2][X].push_back(J[1][0] / surface_J[2][0] * cw);

		surface_normal[0][Y].push_back((J[0][0] - J[0][1]) / surface_J[0][0] * cw);
		surface_normal[1][Y].push_back(J[0][1] / surface_J[1][0] * cw);
		surface_normal[2][Y].push_back(-J[0][0] / surface_J[2][0] * cw);

        return surface_normal;
}
}