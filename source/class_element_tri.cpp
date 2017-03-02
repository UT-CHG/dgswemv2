#include <iostream>
#include <cmath>

#include "class_element_tri.h"

ELEMENT_TRI::ELEMENT_TRI(int ID, int p, int p_geom, double nodal_coordinates_x[],
    double nodal_coordinates_y[], BASIS_TRI* basis, BASIS_GEOM_TRI* basis_geom)
  : ELEMENT(ID,p), p_geom(p_geom), nodal_coordinates_x(nodal_coordinates_x),
    nodal_coordinates_y(nodal_coordinates_y), basis(basis)
{
    if (this->p_geom != 1) {
        if (basis_geom != nullptr) {
            this->basis_geom = basis_geom;
        }
        else {
            printf("\n");
            printf("ELEMENT_TRI - Fatal error!\n");
            printf("NO GEOMETRY BASIS SPECIFIED!\n");
            exit(1);
        }
    }

    this->ComputeGeometry();

    this->ComputeIntegrationFactors();
}

ELEMENT_TRI::~ELEMENT_TRI() {
    delete[] det_J_area;

    delete[] J_inv_t_area[0][0];
    delete[] J_inv_t_area[0][1];
    delete[] J_inv_t_area[1][0];
    delete[] J_inv_t_area[1][1];

    delete[] J_inv_t_area[0];
    delete[] J_inv_t_area[1];

    delete[] J_inv_t_area;

    for (int i = 0; i < 3; i++) {
        delete[] surface_J_edge[i];
        delete[] normal_edge_x[i];
        delete[] normal_edge_y[i];
    }

    delete[] surface_J_edge;
    delete[] normal_edge_x;
    delete[] normal_edge_y;

    for (int i = 0; i < this->basis->GetNumberBasisFunctions(); i++) {
        delete[] area_int_fac_phi[i];
        delete[] area_int_fac_dphidx[i];
        delete[] area_int_fac_dphidy[i];

        for (int j = 0; j < 3; j++) {
            delete[] edge_int_fac_x[j][i];
            delete[] edge_int_fac_y[j][i];
        }
    }

    delete[] area_int_fac_phi;
    delete[] area_int_fac_dphidx;
    delete[] area_int_fac_dphidy;

    for (int i = 0; i < 3; i++) {
        delete[] edge_int_fac_x[i];
        delete[] edge_int_fac_y[i];
    }
    
    delete[] edge_int_fac_x;
    delete[] edge_int_fac_y;
}

void ELEMENT_TRI::ComputeGeometry() {
    if (this->p_geom == 1) {
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
    int number_gp_area = basis->GetIntegrationRuleArea()->GetNumberGP();
    int number_gp_edge = basis->GetIntegrationRuleLine()->GetNumberGP();
    int number_bf = basis->GetNumberBasisFunctions();

    area_int_fac_phi = new double*[number_bf];
    area_int_fac_dphidx = new double*[number_bf];
    area_int_fac_dphidy = new double*[number_bf];
    
    edge_int_fac_x = new double**[3];
    edge_int_fac_y = new double**[3];

    for (int i = 0; i < 3; i++) {
        edge_int_fac_x[i] = new double*[number_bf];
        edge_int_fac_y[i] = new double*[number_bf];
    }

    if (this->p_geom == 1) {
        for (int i = 0; i < number_bf; i++) {
            area_int_fac_phi[i] = new double[number_gp_area];
            area_int_fac_dphidx[i] = new double[number_gp_area];
            area_int_fac_dphidy[i] = new double[number_gp_area];

            for (int j = 0; j < number_gp_area; j++) {
                area_int_fac_phi[i][j] = basis->GetPhiArea()[i][j] *
                    basis->GetIntegrationRuleArea()->GetWeight()[j];

                area_int_fac_dphidx[i][j] = (basis->GetDPhiDZ1Area()[i][j] * J_inv_t_area[0][0][0] +
                    basis->GetDPhiDZ2Area()[i][j] * J_inv_t_area[0][1][0]) *
                    basis->GetIntegrationRuleArea()->GetWeight()[j];

                area_int_fac_dphidy[i][j] = (basis->GetDPhiDZ1Area()[i][j] * J_inv_t_area[1][0][0] +
                    basis->GetDPhiDZ2Area()[i][j] * J_inv_t_area[1][1][0]) *
                    basis->GetIntegrationRuleArea()->GetWeight()[j];
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < number_bf; j++) {
                edge_int_fac_x[i][j] = new double[number_gp_edge];
                edge_int_fac_y[i][j] = new double[number_gp_edge];
                
                for (int k = 0; k < number_gp_edge; k++) {
                    edge_int_fac_x[i][j][k] = basis->GetPhiEdge()[i][j][k] * normal_edge_x[i][0] *
                        basis->GetIntegrationRuleLine()->GetWeight()[k] *
                        surface_J_edge[i][0] / det_J_area[0];
                    
                    edge_int_fac_y[i][j][k] = basis->GetPhiEdge()[i][j][k] * normal_edge_y[i][0] *
                        basis->GetIntegrationRuleLine()->GetWeight()[k] *
                        surface_J_edge[i][0] / det_J_area[0];
                }
            }
        }
    }
    else{
        //Placeholder for cases p_geom > 1
    }
    
    /*
    for (int j = 1; j < 2; j++) {
        printf("%f\t",this->basis->GetIntegrationRuleArea()->GetWeight()[j]);

        printf("%f\t",this->basis->GetPhiArea()[20][j]);
        printf("%f\n",this->area_int_fac_phi[20][j]);

        printf("%f\t",this->basis->GetIntegrationRuleArea()->GetWeight()[j]);

        printf("%f\t",this->J_inv_t_area[0][0][0]);
        printf("%f\t",this->J_inv_t_area[0][1][0]);
        printf("%f\t",this->J_inv_t_area[1][0][0]);
        printf("%f\n",this->J_inv_t_area[1][1][0]);

        printf("%f\t",this->basis->GetDPhiDZ1Area()[20][j]);
        printf("%f\t",this->basis->GetDPhiDZ2Area()[20][j]);

        printf("%f\t",this->area_int_fac_dphidx[20][j]);
        printf("%f\n\n",this->area_int_fac_dphidy[20][j]);
    }
    */
}
