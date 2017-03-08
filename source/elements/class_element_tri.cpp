#include <iostream>
#include <cmath>

#include "../problem/problem_SWE_2D.h"
#include "class_element_tri.h"


ELEMENT_TRI::ELEMENT_TRI(int ID, double nodal_coordinates_x[], double nodal_coordinates_y[],
    BASIS_2D* basis, BASIS_GEOM_2D* basis_geom) : ELEMENT(ID) {

    this->basis = basis;

    int number_bf_geom = 3;

    if (basis_geom != nullptr) {
        this->basis_geom = basis_geom;
        number_bf_geom = this->basis_geom->GetNumberBasisFunctions();
    }

    this->nodal_coordinates_x = new double[number_bf_geom];
    this->nodal_coordinates_y = new double[number_bf_geom];

    for (int i = 0; i < number_bf_geom; i++) {
        this->nodal_coordinates_x[i] = nodal_coordinates_x[i];
        this->nodal_coordinates_y[i] = nodal_coordinates_y[i];
    }

    this->interfaces = new INTERFACE_2D*[3];

    for (int i = 0; i < 3; i++) {
        interfaces[i] = nullptr;
    }

    int number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    int number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    int number_bf = this->basis->GetNumberBasisFunctions();

    this->U = new double*[SIZE_U];
    this->U_area = new double*[SIZE_U];
    this->U_edge = new double**[3];

    for (int i = 0; i < 3; i++) {
        this->U_edge[i] = new double*[SIZE_U];
    }

    for (int i = 0; i < SIZE_U; i++) {
        this->U[i] = new double[number_bf];
        this->U_area[i] = new double[number_gp_area];

        for (int j = 0; j < 3; j++) {
            this->U_edge[j][i] = new double[number_gp_edge];
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
            delete[] this->edge_int_fac_x[j][i];
            delete[] this->edge_int_fac_y[j][i];
        }
    }

    delete[] this->area_int_fac_phi;
    delete[] this->area_int_fac_dphidx;
    delete[] this->area_int_fac_dphidy;

    for (int i = 0; i < 3; i++) {
        delete[] this->edge_int_fac_x[i];
        delete[] this->edge_int_fac_y[i];
    }
    
    delete[] this->edge_int_fac_x;
    delete[] this->edge_int_fac_y;
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
    int number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    int number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    int number_bf = this->basis->GetNumberBasisFunctions();

    this->area_int_fac_phi = new double*[number_bf];
    this->area_int_fac_dphidx = new double*[number_bf];
    this->area_int_fac_dphidy = new double*[number_bf];
    
    this->edge_int_fac_x = new double**[3];
    this->edge_int_fac_y = new double**[3];

    for (int i = 0; i < 3; i++) {
        this->edge_int_fac_x[i] = new double*[number_bf];
        this->edge_int_fac_y[i] = new double*[number_bf];
    }

    if (this->basis_geom == nullptr) {
        for (int i = 0; i < number_bf; i++) {
            this->area_int_fac_phi[i] = new double[number_gp_area];
            this->area_int_fac_dphidx[i] = new double[number_gp_area];
            this->area_int_fac_dphidy[i] = new double[number_gp_area];

            for (int j = 0; j < number_gp_area; j++) {
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
            for (int j = 0; j < number_bf; j++) {
                this->edge_int_fac_x[i][j] = new double[number_gp_edge];
                this->edge_int_fac_y[i][j] = new double[number_gp_edge];
                
                for (int k = 0; k < number_gp_edge; k++) {
                    this->edge_int_fac_x[i][j][k] = this->basis->GetPhiEdge()[i][j][k] * this->normal_edge_x[i][0] *
                        this->basis->GetIntegrationRuleLine()->GetWeight()[k] *
                        this->surface_J_edge[i][0] / this->det_J_area[0];
                    
                    this->edge_int_fac_y[i][j][k] = this->basis->GetPhiEdge()[i][j][k] * this->normal_edge_y[i][0] *
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
    int number_gp_area = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    int number_bf = this->basis->GetNumberBasisFunctions();

    for (int i = 0; i < number_gp_area; i++) {
        this->U_area[u_flag][i] = 0.0;
    }

    for (int i = 0; i < number_bf; i++) {
        for (int j = 0; j < number_gp_area; j++) {
            this->U_area[u_flag][j] += this->U[u_flag][i] * this->basis->GetPhiArea()[i][j];
        }
    }
}

void ELEMENT_TRI::ComputeBoundaryU(int u_flag) {
    int number_gp_edge = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    int number_bf = this->basis->GetNumberBasisFunctions();

    for (int i = 0; i < number_bf; i++) {
        U[u_flag][i] = i;
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < number_gp_edge; j++) {
            this->U_edge[i][u_flag][j] = 0.0;
        }
        
        for (int j = 0; j < number_bf; j++) {
            for (int k = 0; k < number_gp_edge; k++) {
                this->U_edge[i][u_flag][k] += this->U[u_flag][j] * this->basis->GetPhiEdge()[i][j][k];
            }
        }
    }
}

double ELEMENT_TRI::PerformNumericalIntegration(int number_gp, double u[], double int_fac[]) {
    double integral = 0;

    for (int i = 0; i < number_gp; i++) {
        integral += u[i] * int_fac[i];
    }

    return integral;
}