#include "class_element_2D.h"

ELEMENT_2D::ELEMENT_2D(unsigned int ID, unsigned int* neighbor_ID, unsigned char* boundary_type,
	double* nodal_coordinates_x, double* nodal_coordinates_y,
	BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT(ID) {
	this->basis = basis;

	if (basis_geom != nullptr) {
		this->basis_geom = basis_geom;
	}

	//CREATE TRI ELEMENT
	this->Triangle(neighbor_ID, boundary_type, nodal_coordinates_x, nodal_coordinates_y);

	//COMPUTE NUMERICAL INTEGRATION FACTORS
	this->ComputeIntegrationFactors();
}

ELEMENT_2D::~ELEMENT_2D() {
	delete[] this->neighbor_ID;
	delete[] this->boundary_type;

	delete[] this->nodal_coordinates_x;
	delete[] this->nodal_coordinates_y;

	/*for (int i = 0; i < this->number_interfaces; i++) {
	if (this->interface_owner[i]) {
	delete this->interfaces[i];
	}
	}*///at present it deletes all interfaces, however it should keep some for its neighbors
	//delete through mesh class

	delete[] this->interfaces;
	delete[] this->interface_owner;

	for (int i = 0; i < SIZE_U; i++) {
		delete[] this->u[i];
	}

	for (int i = 0; i < SIZE_U_INTERNAL; i++) {
		delete[] this->u_internal[i];
	}

	for (int i = 0; i < this->number_interfaces; i++) {
		for (int j = 0; j < SIZE_U_BOUNDARY; j++) {
			delete[] this->u_boundary[i][j];
		}
		delete[] this->u_boundary[i];
	}

	delete[] this->u;
	delete[] this->u_internal;
	delete[] this->u_boundary;
	
	delete[] this->RHS;

	for (int i = 0; i < this->number_bf; i++) {
		delete[] this->area_int_fac_phi[i];
		delete[] this->area_int_fac_dphidx[i];
		delete[] this->area_int_fac_dphidy[i];

		for (int j = 0; j < this->number_interfaces; j++) {
			delete[] this->edge_int_fac_phi[j][i];
			delete[] this->edge_int_fac_nx[j][i];
			delete[] this->edge_int_fac_ny[j][i];
		}
	}

	delete[] this->area_int_fac_phi;
	delete[] this->area_int_fac_dphidx;
	delete[] this->area_int_fac_dphidy;

	for (int i = 0; i < this->number_interfaces; i++) {
		delete[] this->edge_int_fac_phi[i];
		delete[] this->edge_int_fac_nx[i];
		delete[] this->edge_int_fac_ny[i];
	}

	delete[] this->edge_int_fac_phi;
	delete[] this->edge_int_fac_nx;
	delete[] this->edge_int_fac_ny;
}

void ELEMENT_2D::ComputeIntegrationFactors() {
	this->area_int_fac_phi = new double*[this->number_bf];
	this->area_int_fac_dphidx = new double*[this->number_bf];
	this->area_int_fac_dphidy = new double*[this->number_bf];


	this->edge_int_fac_phi = new double**[this->number_interfaces];
	this->edge_int_fac_nx = new double**[this->number_interfaces];
	this->edge_int_fac_ny = new double**[this->number_interfaces];

	for (int i = 0; i < this->number_interfaces; i++) {
		this->edge_int_fac_phi[i] = new double*[this->number_bf];
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
			this->area_int_fac_phi[i] = new double[this->number_gp_internal];
			this->area_int_fac_dphidx[i] = new double[this->number_gp_internal];
			this->area_int_fac_dphidy[i] = new double[this->number_gp_internal];

			for (int j = 0; j < this->number_gp_internal; j++) {
				this->area_int_fac_phi[i][j] = phi_area[i][j] * w_area[j];

				this->area_int_fac_dphidx[i][j] = (dphi_dz1_area[i][j] * this->J_inv_t_area[0][0][0] +
					dphi_dz2_area[i][j] * this->J_inv_t_area[0][1][0]) * w_area[j];

				this->area_int_fac_dphidy[i][j] = (dphi_dz1_area[i][j] * this->J_inv_t_area[1][0][0] +
					dphi_dz2_area[i][j] * this->J_inv_t_area[1][1][0]) * w_area[j];
			}
		}

		for (int i = 0; i < this->number_interfaces; i++) {
			for (int j = 0; j < this->number_bf; j++) {
				this->edge_int_fac_phi[i][j] = new double[this->number_gp_boundary];
				this->edge_int_fac_nx[i][j] = new double[this->number_gp_boundary];
				this->edge_int_fac_ny[i][j] = new double[this->number_gp_boundary];

				for (int k = 0; k < this->number_gp_boundary; k++) {
					this->edge_int_fac_phi[i][j][k] = phi_edge[i][j][k] *
						w_line[k] * this->surface_J_edge[i][0] / this->det_J_area[0];

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

	//DELETE GEOMETRY DATA RETAIN SURFACE NORMALS
	delete[] this->det_J_area;

	delete[] this->J_inv_t_area[0][0];
	delete[] this->J_inv_t_area[0][1];
	delete[] this->J_inv_t_area[1][0];
	delete[] this->J_inv_t_area[1][1];

	delete[] this->J_inv_t_area[0];
	delete[] this->J_inv_t_area[1];

	delete[] this->J_inv_t_area;

	for (int i = 0; i < this->number_interfaces; i++) {
		delete[] this->surface_J_edge[i];
	}

	delete[] this->surface_J_edge;
}

std::map<unsigned int, INTERFACE*> ELEMENT_2D::CreateInterfaces() {
	std::map<unsigned int, INTERFACE*> internal_interfaces;

	for (int i = 0; i < this->number_interfaces; i++) {
		if (this->interfaces[i] == nullptr) {
			bool straight;
			if (this->basis_geom == nullptr) straight = true;
			else if (this->basis_geom != nullptr) straight = false;

			bool boundary;
			if (this->neighbor_ID[i] == DEFAULT_ID) boundary = true;
			else if (this->neighbor_ID[i] != DEFAULT_ID) boundary = false;

			this->interfaces[i] = new INTERFACE_2D(this->number_gp_boundary, this->u_boundary[i],
				this->normal_edge_x[i], this->normal_edge_y[i], straight, boundary);
			
			this->interface_owner[i] = true;

			if (this->neighbor_ID[i] != DEFAULT_ID) {
				internal_interfaces[this->neighbor_ID[i]] = this->interfaces[i];
			}
		}
	}
	
	//DELETE NO LONGER NEEDED SURFACE NORMALS
	for (int i = 0; i < this->number_interfaces; i++) {
		delete[] this->normal_edge_x[i];
		delete[] this->normal_edge_y[i];
	}

	delete[] this->normal_edge_x;
	delete[] this->normal_edge_y;

	return internal_interfaces;
}

void ELEMENT_2D::AppendInterface(unsigned int neighbor_ID, INTERFACE* interface_ptr) {
	for (int i = 0; i < this->number_interfaces; i++) {
		if (this->neighbor_ID[i] == neighbor_ID) {
			this->interfaces[i] = (INTERFACE_2D*)interface_ptr;

			this->interfaces[i]->SetPointerEX(this->u_boundary[i]);
		}
	}
}

std::vector<std::pair<unsigned char, INTERFACE*>> ELEMENT_2D::GetOwnInterfaces() {
	std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces;

	for (int i = 0; i < this->number_interfaces; i++) {
		if (this->interface_owner[i]) {
			std::pair<unsigned char, INTERFACE*> own_interface;

			own_interface.first = this->boundary_type[i];
			own_interface.second = this->interfaces[i];
			own_interfaces.push_back(own_interface);
		}
	}

	return own_interfaces;
}

void ELEMENT_2D::ComputeInternalU(int u_flag) {
	double** phi_area = this->basis->GetPhiArea();

	for (int i = 0; i < this->number_gp_internal; i++) {
		this->u_internal[u_flag][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_internal; j++) {
			this->u_internal[u_flag][j] += this->u[u_flag][i] * phi_area[i][j];
		}
	}
}

void ELEMENT_2D::ComputeBoundaryU(int u_flag) {
	double*** phi_edge = this->basis->GetPhiEdge();

	for (int k = 0; k < this->number_interfaces; k++) {
		for (int i = 0; i < this->number_gp_boundary; i++) {
			this->u_boundary[k][u_flag][i] = 0.0;
		}

		for (int i = 0; i < this->number_bf; i++) {
			for (int j = 0; j < this->number_gp_boundary; j++) {
				this->u_boundary[k][u_flag][j] += this->u[u_flag][i] * phi_edge[k][i][j];
			}
		}
	}
}

double ELEMENT_2D::IntegrationInternalPhi(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_internal; i++) {
        integral += this->u_internal[u_flag][i] * this->area_int_fac_phi[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationInternalDPhiDX(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_internal; i++) {
        integral += this->u_internal[u_flag][i] * this->area_int_fac_dphidx[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationInternalDPhiDY(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_internal; i++) {
        integral += this->u_internal[u_flag][i] * this->area_int_fac_dphidy[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryPhi(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_interfaces; i++) {
		for (int j = 0; j < this->number_gp_boundary; j++) {
			integral += this->u_boundary[i][u_flag][j] * this->edge_int_fac_phi[i][phi_n][j];
		}
	}

	return integral;
}

double ELEMENT_2D::IntegrationBoundaryNX(int u_flag, int phi_n) {
    double integral = 0;

	for (int i = 0; i < this->number_interfaces; i++) {
		for (int j = 0; j < this->number_gp_boundary; j++) {
			integral += this->u_boundary[i][u_flag][j] * this->edge_int_fac_nx[i][phi_n][j];
		}
	}

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryNY(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_interfaces; i++) {
		for (int j = 0; j < this->number_gp_boundary; j++) {
			integral += this->u_boundary[i][u_flag][j] * this->edge_int_fac_ny[i][phi_n][j];
		}
	}

	return integral;
}

void ELEMENT_2D::SolveLSE(int u_flag) {
	if (this->orthogonal) {
		for (int i = 0; i < this->number_bf; i++) {
			this->u[u_flag][i] = this->m_inv[0][i] * this->RHS[i];
		}
	}
	else if (!(this->orthogonal)) {
		for (int i = 0; i < this->number_bf; i++) {
			this->u[u_flag][i] = 0;
			for (int j = 0; j < this->number_bf; j++) {
				this->u[u_flag][i] += this->m_inv[i][j] * this->RHS[j];
			}
		}
	}
}