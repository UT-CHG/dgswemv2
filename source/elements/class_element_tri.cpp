#include <iostream>
#include <cmath>

#include "../problem/problem_SWE_2D.h"
#include "class_element_tri.h"


ELEMENT_TRI::ELEMENT_TRI(int ID, double nodal_coordinates_x[], double nodal_coordinates_y[],
    BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT(ID) {

    this->basis = basis;

    this->number_bf_geom = 3;

    if (basis_geom != nullptr) {
        this->basis_geom = basis_geom;
        this->number_bf_geom = this->basis_geom->GetNumberBasisFunctions();
    }

    this->nodal_coordinates_x = new double[this->number_bf_geom];
    this->nodal_coordinates_y = new double[this->number_bf_geom];

    for (int i = 0; i < this->number_bf_geom; i++) {
        this->nodal_coordinates_x[i] = nodal_coordinates_x[i];
        this->nodal_coordinates_y[i] = nodal_coordinates_y[i];
    }

    this->interfaces = new INTERFACE_2D*[3];

    for (int i = 0; i < 3; i++) {
        interfaces[i] = nullptr;
    }

    this->number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    this->number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    this->number_bf = this->basis->GetNumberBasisFunctions();

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
    delete[] this->nodal_coordinates_x;
    delete[] this->nodal_coordinates_y;

    for (int i = 0; i < 1; i++) {
        this->interfaces[i]->~INTERFACE_2D();
    }

    delete[] this->interfaces;

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

    for (int i = 0; i < this->basis->GetNumberBasisFunctions(); i++) {
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

        this->J_inv_t_area[0][0][0] =  J[1][1] / (this->det_J_area[0]);
        this->J_inv_t_area[0][1][0] = -J[1][0] / (this->det_J_area[0]);
        this->J_inv_t_area[1][0][0] = -J[0][1] / (this->det_J_area[0]);
        this->J_inv_t_area[1][1][0] =  J[0][0] / (this->det_J_area[0]);

        this->surface_J_edge[0][0] = sqrt(pow((J[0][0] - J[0][1]), 2.0) + pow((J[1][0] - J[1][1]), 2.0));
        this->surface_J_edge[1][0] = sqrt(pow(J[0][1], 2.0) + pow(J[1][1], 2.0));
        this->surface_J_edge[2][0] = sqrt(pow(J[0][0], 2.0) + pow(J[1][0], 2.0));
        
        this->normal_edge_x[0][0] =	-(J[1][0] - J[1][1]) / this->surface_J_edge[0][0];
        this->normal_edge_x[1][0] = -J[1][1] / this->surface_J_edge[1][0];
        this->normal_edge_x[2][0] =  J[1][0] / this->surface_J_edge[2][0];

        this->normal_edge_y[0][0] =  (J[0][0] - J[0][1]) / this->surface_J_edge[0][0];
        this->normal_edge_y[1][0] =  J[0][1] / this->surface_J_edge[1][0];
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
    this->number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    this->number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    this->number_bf = this->basis->GetNumberBasisFunctions();

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
    else{
        //Placeholder for cases p_geom > 1
    }
}

void ELEMENT_TRI::CreateInterfaces() {
    for (int i = 0; i < 3; i++) {
        if (this->interfaces[i] == nullptr) {
            this->interfaces[i] = new INTERFACE_2D(this->U_edge[i]);
            // send this pointer to the corresponding neighbor
        }
    }
}

void ELEMENT_TRI::ComputeInternalU(int u_flag) {
    this->number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    this->number_bf = this->basis->GetNumberBasisFunctions();

    for (int i = 0; i < this->number_gp_area; i++) {
        this->U_area[u_flag][i] = 0.0;
    }

    for (int i = 0; i < this->number_bf; i++) {
        for (int j = 0; j < this->number_gp_area; j++) {
            this->U_area[u_flag][j] += this->U[u_flag][i] * this->basis->GetPhiArea()[i][j];
        }
    }
}

void ELEMENT_TRI::ComputeBoundaryU(int u_flag) {
    this->number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    this->number_bf = this->basis->GetNumberBasisFunctions();

    for (int i = 0; i < this->number_bf; i++) {
        U[u_flag][i] = i;
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < this->number_gp_edge; j++) {
            this->U_edge[i][u_flag][j] = 0.0;
        }
        
        for (int j = 0; j < this->number_bf; j++) {
            for (int k = 0; k < this->number_gp_edge; k++) {
                this->U_edge[i][u_flag][k] += this->U[u_flag][j] * this->basis->GetPhiEdge()[i][j][k];
            }
        }
    }
}

double ELEMENT_TRI::IntegrationInternalPhi(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_area; i++) {
        integral += this->U_area[u_flag][i] * this->area_int_fac_phi[phi_n][i];
    }

    return integral;
}

double ELEMENT_TRI::IntegrationInternalDPhiDX(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_area; i++) {
        integral += this->U_area[u_flag][i] * this->area_int_fac_dphidx[phi_n][i];
    }

    return integral;
}

double ELEMENT_TRI::IntegrationInternalDPhiDY(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_area; i++) {
        integral += this->U_area[u_flag][i] * this->area_int_fac_dphidy[phi_n][i];
    }

    return integral;
}

double ELEMENT_TRI::IntegrationBoundaryNX(int edge_n, int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_edge; i++) {
        integral += this->U_edge[edge_n][u_flag][i] * this->edge_int_fac_nx[edge_n][phi_n][i];
    }

    return integral;
}

double ELEMENT_TRI::IntegrationBoundaryNY(int edge_n, int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_gp_edge; i++) {
		integral += this->U_edge[edge_n][u_flag][i] * this->edge_int_fac_ny[edge_n][phi_n][i];
	}

	return integral;
}