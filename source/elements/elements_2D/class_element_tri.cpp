#include "class_element_tri.h"

ELEMENT_TRI::ELEMENT_TRI(unsigned int ID, unsigned int neighbor_ID[], unsigned char boundary_type[],
	double nodal_coordinates_x[], double nodal_coordinates_y[],
	BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT_2D(ID, basis, basis_geom)
{
	this->number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
	this->number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
	this->number_bf = this->basis->GetNumberBasisFunctions();
	this->number_bf_geom = 3;

	if (basis_geom != nullptr) {
		this->number_bf_geom = this->basis_geom->GetNumberBasisFunctions();
	}

	this->nodal_coordinates_x = new double[this->number_bf_geom];
	this->nodal_coordinates_y = new double[this->number_bf_geom];
	
	this->neighbor_ID = new unsigned int[3];
	this->boundary_type = new unsigned char[3];

	for (int i = 0; i < this->number_bf_geom; i++) {
		this->nodal_coordinates_x[i] = nodal_coordinates_x[i];
		this->nodal_coordinates_y[i] = nodal_coordinates_y[i];
		this->neighbor_ID[i] = neighbor_ID[i];
		this->boundary_type[i] = boundary_type[i];
	}

	this->interfaces = new INTERFACE_2D*[3];
	this->interface_owner = new bool[3];

	for (int i = 0; i < 3; i++) {
		this->interfaces[i] = nullptr;
		this->interface_owner[i] = false;
	}

	this->U = new double*[SIZE_U];
	this->U_area = new double*[SIZE_U];
	this->U_edge = new double**[3];

	for (int i = 0; i < 3; i++) {
		this->U_edge[i] = new double*[SIZE_U];
	}

	for (int i = 0; i < SIZE_U; i++) {
		this->U[i] = new double[this->number_bf];
		this->U_area[i] = new double[this->number_gp_area];

		for (int j = 0; j < 3; j++) {
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

	/*for (int i = 0; i < 3; i++) {
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

		for (int j = 0; j < 3; j++) {
			delete[] this->U_edge[j][i];
		}
	}

	for (int i = 0; i < 3; i++) {
		delete[] this->U_edge[i];
	}

	delete[] this->U;
	delete[] this->U_area;
	delete[] this->U_edge;

	for (int i = 0; i < this->number_bf; i++) {
		delete[] this->area_int_fac_phi[i];
		delete[] this->area_int_fac_dphidx[i];
		delete[] this->area_int_fac_dphidy[i];

		for (int j = 0; j < 3; j++) {
			delete[] this->edge_int_fac_nx[j][i];
			delete[] this->edge_int_fac_ny[j][i];
		}
	}

	delete[] this->area_int_fac_phi;
	delete[] this->area_int_fac_dphidx;
	delete[] this->area_int_fac_dphidy;

	for (int i = 0; i < 3; i++) {
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

		this->surface_J_edge = new double*[3];
		this->normal_edge_x = new double*[3];
		this->normal_edge_y = new double*[3];

		for (int i = 0; i < 3; i++) {
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

	this->edge_int_fac_nx = new double**[3];
	this->edge_int_fac_ny = new double**[3];

	for (int i = 0; i < 3; i++) {
		this->edge_int_fac_nx[i] = new double*[this->number_bf];
		this->edge_int_fac_ny[i] = new double*[this->number_bf];
	}

	if (this->basis_geom == nullptr) {
		for (int i = 0; i < this->number_bf; i++) {
			this->area_int_fac_phi[i] = new double[this->number_gp_area];
			this->area_int_fac_dphidx[i] = new double[this->number_gp_area];
			this->area_int_fac_dphidy[i] = new double[this->number_gp_area];

			for (int j = 0; j < this->number_gp_area; j++) {
				this->area_int_fac_phi[i][j] = this->basis->GetPhiArea()[i][j] *
					this->basis->GetIntegrationRuleArea()->GetWeight()[j];

				this->area_int_fac_dphidx[i][j] = (this->basis->GetDPhiDZ1Area()[i][j] * this->J_inv_t_area[0][0][0] +
					this->basis->GetDPhiDZ2Area()[i][j] * this->J_inv_t_area[0][1][0]) *
					this->basis->GetIntegrationRuleArea()->GetWeight()[j];

				this->area_int_fac_dphidy[i][j] = (this->basis->GetDPhiDZ1Area()[i][j] * this->J_inv_t_area[1][0][0] +
					this->basis->GetDPhiDZ2Area()[i][j] * this->J_inv_t_area[1][1][0]) *
					this->basis->GetIntegrationRuleArea()->GetWeight()[j];
			}
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < this->number_bf; j++) {
				this->edge_int_fac_nx[i][j] = new double[this->number_gp_edge];
				this->edge_int_fac_ny[i][j] = new double[this->number_gp_edge];

				for (int k = 0; k < this->number_gp_edge; k++) {
					this->edge_int_fac_nx[i][j][k] = this->basis->GetPhiEdge()[i][j][k] * this->normal_edge_x[i][0] *
						this->basis->GetIntegrationRuleLine()->GetWeight()[k] *
						this->surface_J_edge[i][0] / this->det_J_area[0];

					this->edge_int_fac_ny[i][j][k] = this->basis->GetPhiEdge()[i][j][k] * this->normal_edge_y[i][0] *
						this->basis->GetIntegrationRuleLine()->GetWeight()[k] *
						this->surface_J_edge[i][0] / this->det_J_area[0];
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

	for (int i = 0; i < 3; i++) {
		delete[] this->surface_J_edge[i];
		delete[] this->normal_edge_x[i];
		delete[] this->normal_edge_y[i];
	}

	delete[] this->surface_J_edge;
	delete[] this->normal_edge_x;
	delete[] this->normal_edge_y;
}

std::map<unsigned int, INTERFACE*> ELEMENT_TRI::CreateInterfaces() {
	std::map<unsigned int, INTERFACE*> internal_interfaces;
	
	for (int i = 0; i < 3; i++) {
		if (this->interfaces[i] == nullptr) {
			this->interfaces[i] = new INTERFACE_2D(this->U_edge[i]);
			this->interface_owner[i] = true;

			if (this->neighbor_ID[i] != DEFAULT_ID) {
				internal_interfaces[this->neighbor_ID[i]] = this->interfaces[i];
			}
		}
	}

	return internal_interfaces;
}

void ELEMENT_TRI::AppendInterface(unsigned int neighbor_ID, INTERFACE* interface_ptr) {
	for (int i = 0; i < 3; i++) {
		if (this->neighbor_ID[i] == neighbor_ID) {
			this->interfaces[i] = (INTERFACE_2D*)interface_ptr;

			this->interfaces[i]->SetPointerEX(this->U_edge[i]);
		}
	}
}

std::vector<std::pair<unsigned char, INTERFACE*>> ELEMENT_TRI::GetOwnInterfaces() {
	std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces;

	for (int i = 0; i < 3; i++) {
		if (this->interface_owner[i]) {
			std::pair<unsigned char, INTERFACE*> own_interface;

			own_interface.first = this->boundary_type[i];
			own_interface.second = this->interfaces[i];
			own_interfaces.push_back(own_interface);
		}
	}

	return own_interfaces;
}