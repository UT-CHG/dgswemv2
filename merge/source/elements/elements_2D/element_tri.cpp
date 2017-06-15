#include "../../class_element.h"
#include "../../class_master_element.h"


void ELEMENT::Triangle() {
	//COMPUTE GEOMETRY AND M_INV
	if (this->basis_geom == nullptr) {
		Array2D<double> J(2);
		J[0].reserve(2);
		J[1].reserve(2);

		J[0].push_back((this->nodal_coordinates[X][1] - this->nodal_coordinates[X][0]) / 2.0);
		J[0].push_back((this->nodal_coordinates[X][2] - this->nodal_coordinates[X][0]) / 2.0);
		J[1].push_back((this->nodal_coordinates[Y][1] - this->nodal_coordinates[Y][0]) / 2.0);
		J[1].push_back((this->nodal_coordinates[Y][2] - this->nodal_coordinates[Y][0]) / 2.0);

		this->det_J_internal.push_back(J[0][0] * J[1][1] - J[0][1] * J[1][0]);

		this->J_inv_internal[0][0].push_back(J[1][1] / (this->det_J_internal[0]));
		this->J_inv_internal[0][1].push_back(-J[0][1] / (this->det_J_internal[0]));
		this->J_inv_internal[1][0].push_back(-J[1][0] / (this->det_J_internal[0]));
		this->J_inv_internal[1][1].push_back(J[0][0] / (this->det_J_internal[0]));

		//printf("%f %f\n%f %f\n\n", this->J_inv_t_area[0][0][0], this->J_inv_t_area[0][1][0],
		//	this->J_inv_t_area[1][0][0], this->J_inv_t_area[1][1][0]);

		this->surface_J_boundary[0].push_back(sqrt(pow((J[0][0] - J[0][1]), 2.0) + pow((J[1][0] - J[1][1]), 2.0)));
		this->surface_J_boundary[1].push_back(sqrt(pow(J[0][1], 2.0) + pow(J[1][1], 2.0)));
		this->surface_J_boundary[2].push_back(sqrt(pow(J[0][0], 2.0) + pow(J[1][0], 2.0)));

		double cw = this->det_J_internal[0] / abs(this->det_J_internal[0]); //CW or CCW

		this->normal_boundary[0][X].push_back((J[1][1] - J[1][0]) / this->surface_J_boundary[0][0] * cw);
		this->normal_boundary[1][X].push_back(-J[1][1] / this->surface_J_boundary[1][0] * cw);
		this->normal_boundary[2][X].push_back(J[1][0] / this->surface_J_boundary[2][0] * cw);

		this->normal_boundary[0][Y].push_back((J[0][0] - J[0][1]) / this->surface_J_boundary[0][0] * cw);
		this->normal_boundary[1][Y].push_back(J[0][1] / this->surface_J_boundary[1][0] * cw);
		this->normal_boundary[2][Y].push_back(-J[0][0] / this->surface_J_boundary[2][0] * cw);
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}


void ELEMENT::InitializeVTKTriangle(std::vector<Point<3>>& points, Array2D<unsigned int>& cells) {
	unsigned int number_pt = points.size();

	double z1;
	double z2;
	double dz = 2.0 / N_DIV;

	if (this->basis_geom == nullptr) {
		for (int i = 0; i <= N_DIV; i++) {
			for (int j = 0; j <= N_DIV - i; j++) {
				points.push_back({0,0,0});

				z1 = -1.0 + dz*j;
				z2 = -1.0 + dz*i;

				points.back()[0] = this->nodal_coordinates[X][0] * (-(z1 + z2) / 2.0) +
					this->nodal_coordinates[X][1] * ((1 + z1) / 2.0) +
					this->nodal_coordinates[X][2] * ((1 + z2) / 2.0);

				points.back()[1] = this->nodal_coordinates[Y][0] * (-(z1 + z2) / 2.0) +
					this->nodal_coordinates[Y][1] * ((1 + z1) / 2.0) +
					this->nodal_coordinates[Y][2] * ((1 + z2) / 2.0);

				points.back()[2] = 0;
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
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


void ELEMENT::WriteCellDataVTKTriangle(std::vector<double>& cell_data, int u_flag) {
	double* u_post = new double[N_DIV*N_DIV];

	for (int i = 0; i < N_DIV*N_DIV; i++) {
		u_post[i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < N_DIV*N_DIV; j++) {
			u_post[j] += this->u[u_flag][i] * this->phi_postprocessor_cell[i][j];
		}
	}

	for (int i = 0; i < N_DIV*N_DIV; i++) {
		cell_data.push_back(u_post[i]);
	}

	delete[] u_post;
}


void ELEMENT::WritePointDataVTKTriangle(std::vector<double>& point_data, int u_flag) {
	double* u_post = new double[(N_DIV + 1)*(N_DIV + 2) / 2];

	for (int i = 0; i < (N_DIV + 1)*(N_DIV + 2) / 2; i++) {
		u_post[i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < (N_DIV + 1)*(N_DIV + 2) / 2; j++) {
			u_post[j] += this->u[u_flag][i] * this->phi_postprocessor_point[i][j];
		}
	}

	for (int i = 0; i < (N_DIV + 1)*(N_DIV + 2) / 2; i++) {
		point_data.push_back(u_post[i]);
	}

	delete[] u_post;
}