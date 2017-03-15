#include "class_element_tri.h"

ELEMENT_TRI::ELEMENT_TRI(unsigned int ID, unsigned int neighbor_ID[], unsigned char boundary_type[],
	double nodal_coordinates_x[], double nodal_coordinates_y[],
	BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT_2D(ID, basis, basis_geom)
{
	this->number_edges = 3;

	this->number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
	this->number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
	this->number_bf = this->basis->GetNumberBasisFunctions();
	this->number_bf_geom = 3;

	if (basis_geom != nullptr) {
		this->number_bf_geom = this->basis_geom->GetNumberBasisFunctions();
	}

	this->nodal_coordinates_x = new double[this->number_bf_geom];
	this->nodal_coordinates_y = new double[this->number_bf_geom];
	
	this->neighbor_ID = new unsigned int[this->number_edges];
	this->boundary_type = new unsigned char[this->number_edges];

	for (int i = 0; i < this->number_bf_geom; i++) {
		this->nodal_coordinates_x[i] = nodal_coordinates_x[i];
		this->nodal_coordinates_y[i] = nodal_coordinates_y[i];
		this->neighbor_ID[i] = neighbor_ID[i];
		this->boundary_type[i] = boundary_type[i];
	}

	this->interfaces = new INTERFACE_2D*[this->number_edges];
	this->interface_owner = new bool[this->number_edges];

	for (int i = 0; i < this->number_edges; i++) {
		this->interfaces[i] = nullptr;
		this->interface_owner[i] = false;
	}

	this->U = new double*[SIZE_U];
	this->U_area = new double*[SIZE_U];
	this->U_edge = new double**[this->number_edges];

	for (int i = 0; i < this->number_edges; i++) {
		this->U_edge[i] = new double*[SIZE_U];
	}

	for (int i = 0; i < SIZE_U; i++) {
		this->U[i] = new double[this->number_bf];
		this->U_area[i] = new double[this->number_gp_area];

		for (int j = 0; j < this->number_edges; j++) {
			this->U_edge[j][i] = new double[this->number_gp_edge];
		}
	}

	this->ComputeGeometry();

	this->ComputeIntegrationFactors();
}

ELEMENT_TRI::~ELEMENT_TRI() {
	delete[] this->neighbor_ID;
	delete[] this->boundary_type;

	delete[] this->nodal_coordinates_x;
	delete[] this->nodal_coordinates_y;

	/*for (int i = 0; i < this->number_edges; i++) {
		if (this->interface_owner[i]) {
			delete this->interfaces[i];
		}
	}*///at present it deletes all interfaces, however it should keep some for its neighbors
	//delete through mesh class

	delete[] this->interfaces;
	delete[] this->interface_owner;

	for (int i = 0; i < SIZE_U; i++) {
		delete[] this->U[i];
		delete[] this->U_area[i];

		for (int j = 0; j < this->number_edges; j++) {
			delete[] this->U_edge[j][i];
		}
	}

	for (int i = 0; i < this->number_edges; i++) {
		delete[] this->U_edge[i];
	}

	delete[] this->U;
	delete[] this->U_area;
	delete[] this->U_edge;

	for (int i = 0; i < this->number_bf; i++) {
		delete[] this->area_int_fac_phi[i];
		delete[] this->area_int_fac_dphidx[i];
		delete[] this->area_int_fac_dphidy[i];

		for (int j = 0; j < this->number_edges; j++) {
			delete[] this->edge_int_fac_nx[j][i];
			delete[] this->edge_int_fac_ny[j][i];
		}
	}

	delete[] this->area_int_fac_phi;
	delete[] this->area_int_fac_dphidx;
	delete[] this->area_int_fac_dphidy;

	for (int i = 0; i < this->number_edges; i++) {
		delete[] this->edge_int_fac_nx[i];
		delete[] this->edge_int_fac_ny[i];
	}

	delete[] this->edge_int_fac_nx;
	delete[] this->edge_int_fac_ny;
}

void ELEMENT_TRI::ComputeGeometry() {
	if (this->basis_geom == nullptr) {
		double** J = new double*[2];
		J[0] = new double[2];
		J[1] = new double[2];

		this->det_J_area = new double[1];

		this->J_inv_t_area = new double**[2];

		this->J_inv_t_area[0] = new double*[2];
		this->J_inv_t_area[1] = new double*[2];

		this->J_inv_t_area[0][1] = new double[1];
		this->J_inv_t_area[0][0] = new double[1];
		this->J_inv_t_area[1][0] = new double[1];
		this->J_inv_t_area[1][1] = new double[1];

		this->surface_J_edge = new double*[this->number_edges];
		this->normal_edge_x = new double*[this->number_edges];
		this->normal_edge_y = new double*[this->number_edges];

		for (int i = 0; i < this->number_edges; i++) {
			this->surface_J_edge[i] = new double[1];
			this->normal_edge_x[i] = new double[1];
			this->normal_edge_y[i] = new double[1];
		}

		J[0][0] = (this->nodal_coordinates_x[1] - this->nodal_coordinates_x[0]) / 2.0;
		J[0][1] = (this->nodal_coordinates_x[2] - this->nodal_coordinates_x[0]) / 2.0;
		J[1][0] = (this->nodal_coordinates_y[1] - this->nodal_coordinates_y[0]) / 2.0;
		J[1][1] = (this->nodal_coordinates_y[2] - this->nodal_coordinates_y[0]) / 2.0;

		this->det_J_area[0] = J[0][0] * J[1][1] - J[0][1] * J[1][0];

		this->J_inv_t_area[0][0][0] = J[1][1] / (this->det_J_area[0]);
		this->J_inv_t_area[0][1][0] = -J[1][0] / (this->det_J_area[0]);
		this->J_inv_t_area[1][0][0] = -J[0][1] / (this->det_J_area[0]);
		this->J_inv_t_area[1][1][0] = J[0][0] / (this->det_J_area[0]);

		this->surface_J_edge[0][0] = sqrt(pow((J[0][0] - J[0][1]), 2.0) + pow((J[1][0] - J[1][1]), 2.0));
		this->surface_J_edge[1][0] = sqrt(pow(J[0][1], 2.0) + pow(J[1][1], 2.0));
		this->surface_J_edge[2][0] = sqrt(pow(J[0][0], 2.0) + pow(J[1][0], 2.0));

		this->normal_edge_x[0][0] = -(J[1][0] - J[1][1]) / this->surface_J_edge[0][0];
		this->normal_edge_x[1][0] = -J[1][1] / this->surface_J_edge[1][0];
		this->normal_edge_x[2][0] = J[1][0] / this->surface_J_edge[2][0];

		this->normal_edge_y[0][0] = (J[0][0] - J[0][1]) / this->surface_J_edge[0][0];
		this->normal_edge_y[1][0] = J[0][1] / this->surface_J_edge[1][0];
		this->normal_edge_y[2][0] = -J[0][0] / this->surface_J_edge[2][0];

		delete[] J[0];
		delete[] J[1];

		delete[] J;
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}

void ELEMENT_TRI::ComputeIntegrationFactors() {
	this->area_int_fac_phi = new double*[this->number_bf];
	this->area_int_fac_dphidx = new double*[this->number_bf];
	this->area_int_fac_dphidy = new double*[this->number_bf];

	this->edge_int_fac_nx = new double**[this->number_edges];
	this->edge_int_fac_ny = new double**[this->number_edges];

	for (int i = 0; i < this->number_edges; i++) {
		this->edge_int_fac_nx[i] = new double*[this->number_bf];
		this->edge_int_fac_ny[i] = new double*[this->number_bf];
	}

	double** phi_area = this->basis->GetPhiArea();
	double** dphi_dz1_area = this->basis->GetDPhiDZ1Area();
	double** dphi_dz2_area = this->basis->GetDPhiDZ2Area();
	double*** phi_edge = this->basis->GetPhiEdge();

	double* w_area = this->basis->GetIntegrationRuleArea()->GetWeight();
	double* w_line = this->basis->GetIntegrationRuleLine()->GetWeight();

	if (this->basis_geom == nullptr) {
		for (int i = 0; i < this->number_bf; i++) {
			this->area_int_fac_phi[i] = new double[this->number_gp_area];
			this->area_int_fac_dphidx[i] = new double[this->number_gp_area];
			this->area_int_fac_dphidy[i] = new double[this->number_gp_area];

			for (int j = 0; j < this->number_gp_area; j++) {
				this->area_int_fac_phi[i][j] = phi_area[i][j] * w_area[j];

				this->area_int_fac_dphidx[i][j] = (dphi_dz1_area[i][j] * this->J_inv_t_area[0][0][0] +
					dphi_dz2_area[i][j] * this->J_inv_t_area[0][1][0]) * w_area[j];

				this->area_int_fac_dphidy[i][j] = (dphi_dz1_area[i][j] * this->J_inv_t_area[1][0][0] +
					dphi_dz2_area[i][j] * this->J_inv_t_area[1][1][0]) * w_area[j];
			}
		}

		for (int i = 0; i < this->number_edges; i++) {
			for (int j = 0; j < this->number_bf; j++) {
				this->edge_int_fac_nx[i][j] = new double[this->number_gp_edge];
				this->edge_int_fac_ny[i][j] = new double[this->number_gp_edge];

				for (int k = 0; k < this->number_gp_edge; k++) {
					this->edge_int_fac_nx[i][j][k] = phi_edge[i][j][k] * this->normal_edge_x[i][0] *
						w_line[k] * this->surface_J_edge[i][0] / this->det_J_area[0];

					this->edge_int_fac_ny[i][j][k] = phi_edge[i][j][k] * this->normal_edge_y[i][0] *
						w_line[k] * this->surface_J_edge[i][0] / this->det_J_area[0];
				}
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
	}

	// Delete geometry data
	delete[] this->det_J_area;

	delete[] this->J_inv_t_area[0][0];
	delete[] this->J_inv_t_area[0][1];
	delete[] this->J_inv_t_area[1][0];
	delete[] this->J_inv_t_area[1][1];

	delete[] this->J_inv_t_area[0];
	delete[] this->J_inv_t_area[1];

	delete[] this->J_inv_t_area;

	for (int i = 0; i < this->number_edges; i++) {
		delete[] this->surface_J_edge[i];
		delete[] this->normal_edge_x[i];
		delete[] this->normal_edge_y[i];
	}

	delete[] this->surface_J_edge;
	delete[] this->normal_edge_x;
	delete[] this->normal_edge_y;
}