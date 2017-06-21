#include "../shapes_2D.h"

namespace Shape {
	std::vector<double> StraightTriangle::get_J_det(const std::vector<Point<2>>& nodal_coordinates) {
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

	Array3D<double> StraightTriangle::get_J_inv(const std::vector<Point<2>>& nodal_coordinates) {
		Array3D<double> J_inv(2);
		J_inv[0].resize(2);
		J_inv[1].resize(2);

		Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((nodal_coordinates[1][X] - nodal_coordinates[0][X]) / 2.0);
		J[0].push_back((nodal_coordinates[2][X] - nodal_coordinates[0][X]) / 2.0);
		J[1].push_back((nodal_coordinates[1][Y] - nodal_coordinates[0][Y]) / 2.0);
		J[1].push_back((nodal_coordinates[2][Y] - nodal_coordinates[0][Y]) / 2.0);

		double det_J = this->get_J_det(nodal_coordinates)[0];

		J_inv[0][0].push_back(J[1][1] / det_J);
		J_inv[0][1].push_back(-J[0][1] / det_J);
		J_inv[1][0].push_back(-J[1][0] / det_J);
		J_inv[1][1].push_back(J[0][0] / det_J);

		return J_inv;
	}

	std::vector<double> StraightTriangle::get_surface_J(int n_bound, const std::vector<Point<2>>& nodal_coordinates) {
		std::vector<double> surface_J;

		int pt_begin, pt_end;

		if (n_bound == 0) {
			pt_begin = 1;
			pt_end = 2;
		}
		else if (n_bound == 1) {
			pt_begin = 2;
			pt_end = 0;
		}
		else if (n_bound == 2) {
			pt_begin = 0;
			pt_end = 1;
		}

		surface_J.push_back(sqrt(pow(nodal_coordinates[pt_end][X] - nodal_coordinates[pt_begin][X], 2.0) +
			pow(nodal_coordinates[pt_end][Y] - nodal_coordinates[pt_begin][Y], 2.0)) / 2.0); //half length for straight edge

		return surface_J;
	}

	Array2D<double> StraightTriangle::get_surface_normal(int n_bound, const std::vector<Point<2>>& nodal_coordinates) {
		Array2D<double> surface_normal(1);

		Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((nodal_coordinates[1][X] - nodal_coordinates[0][X]) / 2.0);
		J[0].push_back((nodal_coordinates[2][X] - nodal_coordinates[0][X]) / 2.0);
		J[1].push_back((nodal_coordinates[1][Y] - nodal_coordinates[0][Y]) / 2.0);
		J[1].push_back((nodal_coordinates[2][Y] - nodal_coordinates[0][Y]) / 2.0);

		double det_J = (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
		double cw = det_J / std::abs(det_J); //CW or CCW

		int pt_begin, pt_end;

		if (n_bound == 0) {
			pt_begin = 1;
			pt_end = 2;
		}
		else if (n_bound == 1) {
			pt_begin = 2;
			pt_end = 0;
		}
		else if (n_bound == 2) {
			pt_begin = 0;
			pt_end = 1;
		}

		double length = sqrt(pow(nodal_coordinates[pt_end][X] - nodal_coordinates[pt_begin][X], 2.0) +
			pow(nodal_coordinates[pt_end][Y] - nodal_coordinates[pt_begin][Y], 2.0));

		surface_normal[0].push_back(cw * (nodal_coordinates[pt_end][Y] - nodal_coordinates[pt_begin][Y]) / length);
		surface_normal[0].push_back(-cw * (nodal_coordinates[pt_end][X] - nodal_coordinates[pt_begin][X]) / length);

		return surface_normal;
	}

	void StraightTriangle::get_VTK(std::vector<Point<3>>& points, Array2D<unsigned int>& cells, const std::vector<Point<2>>& nodal_coordinates) {
		unsigned int number_pt = points.size();

		double z1;
		double z2;
		double dz = 2.0 / N_DIV;

		for (int i = 0; i <= N_DIV; i++) {
			for (int j = 0; j <= N_DIV - i; j++) {
				points.push_back({ 0,0,0 });

				z1 = -1.0 + dz*j;
				z2 = -1.0 + dz*i;

				points.back()[0] = nodal_coordinates[0][X] * (-(z1 + z2) / 2.0) +
					nodal_coordinates[1][X] * ((1 + z1) / 2.0) +
					nodal_coordinates[2][X] * ((1 + z2) / 2.0);

				points.back()[1] = nodal_coordinates[0][Y] * (-(z1 + z2) / 2.0) +
					nodal_coordinates[1][Y] * ((1 + z1) / 2.0) +
					nodal_coordinates[2][Y] * ((1 + z2) / 2.0);

				points.back()[2] = 0;
			}
		}

		unsigned int pt_ID;

		for (int i = 0; i < N_DIV; i++) {
			for (int j = 0; j < N_DIV - i; j++) {
				cells.push_back(std::vector<unsigned int>(4));

				pt_ID = number_pt + (N_DIV + 1)*(N_DIV + 2) / 2 - (N_DIV - i + 1)*(N_DIV - i + 2) / 2 + j;

				cells.back()[0] = 5;
				cells.back()[1] = pt_ID;
				cells.back()[2] = pt_ID + 1;
				cells.back()[3] = pt_ID + (N_DIV + 1 - i);
			}
		}

		for (int i = 1; i < N_DIV; i++) {
			for (int j = 0; j < N_DIV - i; j++) {
				cells.push_back(std::vector<unsigned int>(4));

				pt_ID = number_pt + (N_DIV + 1)*(N_DIV + 2) / 2 - (N_DIV - i + 1)*(N_DIV - i + 2) / 2 + j;

				cells.back()[0] = 5;
				cells.back()[1] = pt_ID;
				cells.back()[2] = pt_ID + 1;
				cells.back()[3] = pt_ID - (N_DIV + 2 - i) + 1;
			}
		}
	}
}