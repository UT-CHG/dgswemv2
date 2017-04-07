#include <cmath>

#include "../class_element_2D.h"

void ELEMENT_2D::Triangle(unsigned int* neighbor_ID, unsigned char* boundary_type,
    double* nodal_coordinates_x, double* nodal_coordinates_y)
{
    this->number_interfaces = 3;

    this->number_gp_internal = this->basis->GetIntegrationRuleArea()->GetNumberGP();
    this->number_gp_boundary = this->basis->GetIntegrationRuleLine()->GetNumberGP();
    this->number_bf = this->basis->GetNumberBasisFunctions();
    this->number_bf_geom = 3;

    if (basis_geom != nullptr) {
        this->number_bf_geom = this->basis_geom->GetNumberBasisFunctions();
    }

    this->nodal_coordinates_x = new double[this->number_bf_geom];
    this->nodal_coordinates_y = new double[this->number_bf_geom];
    
    for (int i = 0; i < this->number_bf_geom; i++) {
        this->nodal_coordinates_x[i] = nodal_coordinates_x[i];
        this->nodal_coordinates_y[i] = nodal_coordinates_y[i];
    }

    this->neighbor_ID = new unsigned int[this->number_interfaces];
    this->boundary_type = new unsigned char[this->number_interfaces];
    this->interfaces = new INTERFACE_2D*[this->number_interfaces];
    this->interface_owner = new bool[this->number_interfaces];

    for (int i = 0; i < this->number_interfaces; i++) {
        this->neighbor_ID[i] = neighbor_ID[i];
        this->boundary_type[i] = boundary_type[i];
        this->interfaces[i] = nullptr;
        this->interface_owner[i] = false;
    }

    //INITIALIZE ARRAYS TO STORE Us and RHS
    this->RHS = new double[this->number_bf];

    this->u = new double*[SIZE_U];
    this->u_internal = new double*[SIZE_U_INTERNAL];
    this->u_boundary = new double**[this->number_interfaces];

    for (int i = 0; i < SIZE_U; i++) {
        this->u[i] = new double[this->number_bf];
    }

    for (int i = 0; i < SIZE_U_INTERNAL; i++) {
        this->u_internal[i] = new double[this->number_gp_internal];
    }

    for (int i = 0; i < this->number_interfaces; i++) {
        this->u_boundary[i] = new double*[SIZE_U_BOUNDARY];
        for (int j = 0; j < SIZE_U_BOUNDARY; j++) {
            this->u_boundary[i][j] = new double[this->number_gp_boundary];
        }
    }

    //COMPUTE GEOMETRY AND M_INV
    if (this->basis_geom == nullptr) {
        this->orthogonal = this->basis->GetOrthogonal();
        this->m_inv = this->basis->GetMInv();

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

        this->surface_J_edge = new double*[this->number_interfaces];
        this->normal_edge_x = new double*[this->number_interfaces];
        this->normal_edge_y = new double*[this->number_interfaces];

        for (int i = 0; i < this->number_interfaces; i++) {
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

        //printf("%f %f\n%f %f\n\n", this->J_inv_t_area[0][0][0], this->J_inv_t_area[0][1][0],
        //	this->J_inv_t_area[1][0][0], this->J_inv_t_area[1][1][0]);

        this->surface_J_edge[0][0] = sqrt(pow((J[0][0] - J[0][1]), 2.0) + pow((J[1][0] - J[1][1]), 2.0));
        this->surface_J_edge[1][0] = sqrt(pow(J[0][1], 2.0) + pow(J[1][1], 2.0));
        this->surface_J_edge[2][0] = sqrt(pow(J[0][0], 2.0) + pow(J[1][0], 2.0));

        double cw = this->det_J_area[0] / abs(this->det_J_area[0]); //CW or CCW

        this->normal_edge_x[0][0] = (J[1][1] - J[1][0]) / this->surface_J_edge[0][0] * cw;
        this->normal_edge_x[1][0] = -J[1][1] / this->surface_J_edge[1][0] * cw;
        this->normal_edge_x[2][0] = J[1][0] / this->surface_J_edge[2][0] *cw;

        this->normal_edge_y[0][0] = (J[0][0] - J[0][1]) / this->surface_J_edge[0][0] * cw;
        this->normal_edge_y[1][0] = J[0][1] / this->surface_J_edge[1][0] * cw;
        this->normal_edge_y[2][0] = -J[0][0] / this->surface_J_edge[2][0] * cw;

        delete[] J[0];
        delete[] J[1];

        delete[] J;
    }
    else {
        //Placeholder for cases p_geom > 1
    }
}

void ELEMENT_2D::InitializeVTKTriangle(std::vector<double*>& points, std::vector<unsigned int*>& cells) {
    unsigned int number_pt = points.size();
    
    if (this->basis_geom == nullptr) {
        double z1;
        double z2;
        double dz = 2.0 / N_DIV;

        for (int i = 0; i <= N_DIV; i++) {
            for (int j = 0; j <= N_DIV - i; j++) {
                points.push_back(new double[3]);

                z1 = -1.0 + dz*j;
                z2 = -1.0 + dz*i;

                points.back()[0] = this->nodal_coordinates_x[0] * (-(z1 + z2) / 2.0) +
                    this->nodal_coordinates_x[1] * ((1 + z1) / 2.0) +
                    this->nodal_coordinates_x[2] * ((1 + z2) / 2.0);

                points.back()[1] = this->nodal_coordinates_y[0] * (-(z1 + z2) / 2.0) +
                    this->nodal_coordinates_y[1] * ((1 + z1) / 2.0) +
                    this->nodal_coordinates_y[2] * ((1 + z2) / 2.0);

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