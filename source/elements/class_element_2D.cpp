#include "class_element_2D.h"

ELEMENT_2D::ELEMENT_2D(unsigned int ID, BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT(ID) {
	this->basis = basis;

	if (basis_geom != nullptr) {
		this->basis_geom = basis_geom;
	}
}

void ELEMENT_2D::ComputeInternalU(int u_flag) {
    for (int i = 0; i < this->number_gp_area; i++) {
        this->U_area[u_flag][i] = 0.0;
    }

    for (int i = 0; i < this->number_bf; i++) {
        for (int j = 0; j < this->number_gp_area; j++) {
            this->U_area[u_flag][j] += this->U[u_flag][i] * this->basis->GetPhiArea()[i][j];
        }
    }
}

void ELEMENT_2D::ComputeBoundaryU(int edge_n, int u_flag) {
	for (int i = 0; i < this->number_gp_edge; i++) {
		this->U_edge[edge_n][u_flag][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_edge; j++) {
			this->U_edge[edge_n][u_flag][j] += this->U[u_flag][i] * this->basis->GetPhiEdge()[edge_n][i][j];
		}
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

double ELEMENT_2D::IntegrationBoundaryNX(int edge_n, int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_edge; i++) {
        integral += this->U_edge[edge_n][u_flag][i] * this->edge_int_fac_nx[edge_n][phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryNY(int edge_n, int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_gp_edge; i++) {
		integral += this->U_edge[edge_n][u_flag][i] * this->edge_int_fac_ny[edge_n][phi_n][i];
	}

	return integral;
}