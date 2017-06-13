#include "../../class_element.h"
#include "../../class_master_element.h"

//template<int element_type = TRIANGLE, class Basis = Dubiner_2D,
//	class Integration_int = Dunavant_2D, class Integration_bound = GaussLegendre_1D >
void MasterElement::MasterTriangle(int p) {
	std::pair<std::vector<double>, std::vector<Point<1>>> integration_rule_boundary = this->integration_boundary.get_rule(2 * p);
	std::pair<std::vector<double>, std::vector<Point<2>>> integration_rule_internal = this->integration_internal.get_rule(2 * p);
	
	this->m_inv = this->basis.get_m_inv(p);
	
	this->phi_internal = this->basis.get_phi(p, integration_rule_internal.second);
	this->basis.basis_test(p, this->phi_internal, integration_rule_internal); //TEST BASIS

	this->dphi_internal = this->basis.get_dphi(p, integration_rule_internal.second);
	
	this->phi_boundary.reserve(3);
	std::vector<Point<2>> z_boundary(integration_rule_boundary.first.size());
	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			for (size_t j = 0; j < z_boundary.size(); j++) {
				z_boundary[j][Z1] = -integration_rule_boundary.second[j][Z1];
				z_boundary[j][Z2] = integration_rule_boundary.second[j][Z1];
			}
		}
		else if (i == 1) {
			for (size_t j = 0; j < z_boundary.size(); j++) {
				z_boundary[j][Z1] = -1;
				z_boundary[j][Z2] = -integration_rule_boundary.second[j][Z1];
			}
		}
		else if (i == 2) {
			for (size_t j = 0; j < z_boundary.size(); j++) {
				z_boundary[j][Z1] = integration_rule_boundary.second[j][Z1];
				z_boundary[j][Z2] = -1;
			}
		}
		this->phi_boundary.push_back(basis.get_phi(p, z_boundary));
	}

	this->internal_int_fac_phi = this->phi_internal;
	for (size_t i = 0; i < this->internal_int_fac_phi.size(); i++) { //iterate through basis functions
		for (size_t j = 0; j < internal_int_fac_phi[i].size(); j++) { //iterate through internal GPs
			this->internal_int_fac_phi[i][j] *= integration_rule_internal.first[j]; //apply weight
		}
	}

	this->internal_int_fac_dphi = this->dphi_internal;
	for (size_t i = 0; i < internal_int_fac_dphi.size(); i++) { //iterate through basis functions
		for (size_t j = 0; j < internal_int_fac_dphi[i].size(); j++) { //iterate through differentiation directions
			for (size_t k = 0; k < internal_int_fac_dphi[i][j].size(); k++) { //iterate through internal GPs
				this->internal_int_fac_dphi[i][j][k] *= integration_rule_internal.first[k]; //apply weight
			}
		}
	}

	this->boundary_int_fac_phi = this->phi_boundary;
	for (size_t i = 0; i < boundary_int_fac_phi.size(); i++) { //iterate thorough boundaries
		for (size_t j = 0; j < boundary_int_fac_phi[i].size(); j++) { //iterate through basis funtions
			for (size_t k = 0; k < boundary_int_fac_phi[i][j].size(); k++) { //iterate through boundary GPs
				this->boundary_int_fac_phi[i][j][k] *= integration_rule_boundary.first[k]; //apply weight
			}
		}
	}

	double dz = 2.0 / N_DIV;

	std::vector<Point<2>> z_postprocessor_cell(N_DIV*N_DIV);
	int n_pt = 0;
	for (int i = 0; i < N_DIV; i++) {
		for (int j = 0; j < N_DIV - i; j++) {
			z_postprocessor_cell[n_pt][Z1] = -1.0 + dz*j + dz / 3.0; //CENTROID
			z_postprocessor_cell[n_pt][Z2] = -1.0 + dz*i + dz / 3.0;
			n_pt++;
		}
	}
	for (int i = 1; i < N_DIV; i++) {
		for (int j = 0; j < N_DIV - i; j++) {
			z_postprocessor_cell[n_pt][Z1] = -1.0 + dz*j + 2 * dz / 3.0; //CENTROID
			z_postprocessor_cell[n_pt][Z2] = -1.0 + dz*i - dz / 3.0;
			n_pt++;
		}
	}

	this->phi_postprocessor_cell = this->basis.get_phi(p, z_postprocessor_cell);

	std::vector<Point<2>> z_postprocessor_point((N_DIV + 1)*(N_DIV + 2) / 2);
	n_pt = 0;
	for (int i = 0; i < N_DIV; i++) {
		for (int j = 0; j <= N_DIV - i; j++) {
			z_postprocessor_point[n_pt][Z1] = -1.0 + dz*j;
			z_postprocessor_point[n_pt][Z2] = -1.0 + dz*i;
			n_pt++;
		}
	}

	this->phi_postprocessor_point = this->basis.get_phi(p, z_postprocessor_point);
}

void ELEMENT::Triangle() {
	this->dimension = 2;
	this->number_boundaries = 3;

	this->number_bf = this->master.phi_internal.size();
	this->number_gp_internal = this->master.phi_internal[0].size();
	this->number_gp_boundary = this->master.phi_boundary[0][0].size();

	if (basis_geom != nullptr) {
		this->number_bf_geom = 3;
	}
	else {
		this->number_bf_geom = 3;
	}

	this->allocate_memory();

	for (int i = 0; i < this->number_boundaries; i++) {
		this->interfaces[i] = nullptr;
		this->interface_owner[i] = false;
	}

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

void ELEMENT::InitializeVTKTriangle(std::vector<double*>& points, std::vector<unsigned int*>& cells) {
	unsigned int number_pt = points.size();

	double z1;
	double z2;
	double dz = 2.0 / N_DIV;

	if (this->basis_geom == nullptr) {
		for (int i = 0; i <= N_DIV; i++) {
			for (int j = 0; j <= N_DIV - i; j++) {
				points.push_back(new double[3]);

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
			cells.push_back(new unsigned int[4]);

			pt_ID = number_pt + (N_DIV + 1)*(N_DIV + 2) / 2 - (N_DIV - i + 1)*(N_DIV - i + 2) / 2 + j;

			cells.back()[0] = 5;
			cells.back()[1] = pt_ID;
			cells.back()[2] = pt_ID + 1;
			cells.back()[3] = pt_ID + (N_DIV + 1 - i);
		}
	}

	for (int i = 1; i < N_DIV; i++) {
		for (int j = 0; j < N_DIV - i; j++) {
			cells.push_back(new unsigned int[4]);

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