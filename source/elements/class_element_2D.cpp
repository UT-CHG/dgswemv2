#include "class_element_2D.h"

ELEMENT_2D::ELEMENT_2D(unsigned int ID, BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT(ID) {
	this->basis = basis;

	if (basis_geom != nullptr) {
		this->basis_geom = basis_geom;
	}
}

std::map<unsigned int, INTERFACE*> ELEMENT_2D::CreateInterfaces() {
	std::map<unsigned int, INTERFACE*> internal_interfaces;

	for (int i = 0; i < this->number_edges; i++) {
		if (this->interfaces[i] == nullptr) {
			this->interfaces[i] = new INTERFACE_2D(this->number_gp_edge, this->U_edge[i]);
			this->interface_owner[i] = true;

			if (this->neighbor_ID[i] != DEFAULT_ID) {
				internal_interfaces[this->neighbor_ID[i]] = this->interfaces[i];
			}
		}
	}

	return internal_interfaces;
}

void ELEMENT_2D::AppendInterface(unsigned int neighbor_ID, INTERFACE* interface_ptr) {
	for (int i = 0; i < this->number_edges; i++) {
		if (this->neighbor_ID[i] == neighbor_ID) {
			this->interfaces[i] = (INTERFACE_2D*)interface_ptr;

			this->interfaces[i]->SetPointerEX(this->U_edge[i]);
		}
	}
}

std::vector<std::pair<unsigned char, INTERFACE*>> ELEMENT_2D::GetOwnInterfaces() {
	std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces;

	for (int i = 0; i < this->number_edges; i++) {
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

	for (int i = 0; i < this->number_gp_area; i++) {
		this->U_area[u_flag][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_area; j++) {
			this->U_area[u_flag][j] += this->U[u_flag][i] * phi_area[i][j];
		}
	}
}

void ELEMENT_2D::ComputeBoundaryU(int u_flag) {
	double*** phi_edge = this->basis->GetPhiEdge();

	for (int k = 0; k < this->number_edges; k++) {
		for (int i = 0; i < this->number_gp_edge; i++) {
			this->U_edge[k][u_flag][i] = 0.0;
		}

		for (int i = 0; i < this->number_bf; i++) {
			for (int j = 0; j < this->number_gp_edge; j++) {
				this->U_edge[k][u_flag][j] += this->U[u_flag][i] * phi_edge[k][i][j];
			}
		}
	}
}

void ELEMENT_2D::ComputeF() {
	problem_compute_f(this->number_gp_area, this->U_area);

	for (int i = 0; i < this->number_edges; i++) {
		problem_compute_f(this->number_gp_edge, this->U_edge[i]);
	}
}

double ELEMENT_2D::IntegrationInternalPhi(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_area; i++) {
        integral += this->U_area[u_flag][i] * this->area_int_fac_phi[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationInternalDPhiDX(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_area; i++) {
        integral += this->U_area[u_flag][i] * this->area_int_fac_dphidx[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationInternalDPhiDY(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_area; i++) {
        integral += this->U_area[u_flag][i] * this->area_int_fac_dphidy[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryNX(int u_flag, int phi_n) {
    double integral = 0;

	for (int i = 0; i < this->number_edges; i++) {
		for (int j = 0; j < this->number_gp_edge; j++) {
			integral += this->U_edge[i][u_flag][j] * this->edge_int_fac_nx[i][phi_n][j];
		}
	}

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryNY(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_edges; i++) {
		for (int j = 0; j < this->number_gp_edge; j++) {
			integral += this->U_edge[i][u_flag][j] * this->edge_int_fac_ny[i][phi_n][j];
		}
	}

	return integral;
}