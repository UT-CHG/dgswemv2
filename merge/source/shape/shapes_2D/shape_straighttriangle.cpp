#include "../shapes_2D.hpp"

namespace Shape {
	std::vector<double> StraightTriangle::GetJdet(const std::vector<Point<2>>& pts) {
		std::vector<double> J_det;

		Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
		J[0].push_back((this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
		J[1].push_back((this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);
		J[1].push_back((this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);

		J_det.push_back(J[0][0] * J[1][1] - J[0][1] * J[1][0]);

		return J_det;
	}

	Array3D<double> StraightTriangle::GetJinv(const std::vector<Point<2>>& pts) {
		Array3D<double> J_inv(2);
		J_inv[0].resize(2);
		J_inv[1].resize(2);

		Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
		J[0].push_back((this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
		J[1].push_back((this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);
		J[1].push_back((this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);

		double det_J = (J[0][0] * J[1][1] - J[0][1] * J[1][0]);

		J_inv[0][0].push_back(J[1][1] / det_J);
		J_inv[0][1].push_back(-J[0][1] / det_J);
		J_inv[1][0].push_back(-J[1][0] / det_J);
		J_inv[1][1].push_back(J[0][0] / det_J);

		return J_inv;
	}

	std::vector<double> StraightTriangle::GetSurfaceJ(uint n_bound, const std::vector<Point<2>>& pts) {
		std::vector<double> surface_J;

		uint pt_begin, pt_end;

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

		surface_J.push_back(sqrt(pow(this->nodal_coordinates[pt_end][GlobalCoord::x] - this->nodal_coordinates[pt_begin][GlobalCoord::x], 2.0) +
			pow(this->nodal_coordinates[pt_end][GlobalCoord::y] - this->nodal_coordinates[pt_begin][GlobalCoord::y], 2.0)) / 2.0); //half length for straight edge

		return surface_J;
	}

	Array2D<double> StraightTriangle::GetSurfaceNormal(uint n_bound, const std::vector<Point<2>>& pts) {
		Array2D<double> surface_normal(1);

		Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
		J[0].push_back((this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
		J[1].push_back((this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);
		J[1].push_back((this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);

		double det_J = (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
		double cw = det_J / std::abs(det_J); //CW or CCW

		uint pt_begin, pt_end;

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

		double length = sqrt(pow(this->nodal_coordinates[pt_end][GlobalCoord::x] - this->nodal_coordinates[pt_begin][GlobalCoord::x], 2.0) +
			pow(this->nodal_coordinates[pt_end][GlobalCoord::y] - this->nodal_coordinates[pt_begin][GlobalCoord::y], 2.0));

		surface_normal[0].push_back(cw * (this->nodal_coordinates[pt_end][GlobalCoord::y] - this->nodal_coordinates[pt_begin][GlobalCoord::y]) / length);
		surface_normal[0].push_back(-cw * (this->nodal_coordinates[pt_end][GlobalCoord::x] - this->nodal_coordinates[pt_begin][GlobalCoord::x]) / length);

		return surface_normal;
	}

	void StraightTriangle::GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) {
		uint number_pt = points.size();

		double z1;
		double z2;
		double dz = 2.0 / N_DIV;

		for (uint i = 0; i <= N_DIV; i++) {
			for (uint j = 0; j <= N_DIV - i; j++) {
				points.push_back({ 0,0,0 });

				z1 = -1.0 + dz*j;
				z2 = -1.0 + dz*i;

				points.back()[0] = this->nodal_coordinates[0][GlobalCoord::x] * (-(z1 + z2) / 2.0) +
					this->nodal_coordinates[1][GlobalCoord::x] * ((1 + z1) / 2.0) +
					this->nodal_coordinates[2][GlobalCoord::x] * ((1 + z2) / 2.0);

				points.back()[1] = this->nodal_coordinates[0][GlobalCoord::y] * (-(z1 + z2) / 2.0) +
					this->nodal_coordinates[1][GlobalCoord::y] * ((1 + z1) / 2.0) +
					this->nodal_coordinates[2][GlobalCoord::y] * ((1 + z2) / 2.0);

				points.back()[2] = 0;
			}
		}

		uint pt_ID;

		for (uint i = 0; i < N_DIV; i++) {
			for (uint j = 0; j < N_DIV - i; j++) {
				cells.push_back(std::vector<uint>(4));

				pt_ID = number_pt + (N_DIV + 1)*(N_DIV + 2) / 2 - (N_DIV - i + 1)*(N_DIV - i + 2) / 2 + j;

				cells.back()[0] = 5;
				cells.back()[1] = pt_ID;
				cells.back()[2] = pt_ID + 1;
				cells.back()[3] = pt_ID + (N_DIV + 1 - i);
			}
		}

		for (uint i = 1; i < N_DIV; i++) {
			for (uint j = 0; j < N_DIV - i; j++) {
				cells.push_back(std::vector<uint>(4));

				pt_ID = number_pt + (N_DIV + 1)*(N_DIV + 2) / 2 - (N_DIV - i + 1)*(N_DIV - i + 2) / 2 + j;

				cells.back()[0] = 5;
				cells.back()[1] = pt_ID;
				cells.back()[2] = pt_ID + 1;
				cells.back()[3] = pt_ID - (N_DIV + 2 - i) + 1;
			}
		}
	}

	std::vector<double> StraightTriangle::InterpolateNodalValues(const std::vector<double>& nodal_values, const std::vector<Point<2>>& pts) {
		std::vector<double> interpolation;
		
		interpolation.reserve(pts.size());

		for (uint pt = 0; pt < pts.size(); pt++) { 
			interpolation.push_back(0);
			
			interpolation[pt] = 
				-(pts[pt][LocalCoordTri::z1] + pts[pt][LocalCoordTri::z2])/2 * nodal_values[0] //N1
				+(1 + pts[pt][LocalCoordTri::z1])/2 * nodal_values[1]						   //N2
				+(1 + pts[pt][LocalCoordTri::z2])/2 * nodal_values[2];						   //N3
		}

		return interpolation; 
	}
}