#include <cmath>

#include "../class_element.h"

void ELEMENT::Triangle(unsigned int* neighbor_ID, unsigned char* boundary_type,
    double** nodal_coordinates)
{
	this->dimension = 2;
	this->element_type = TRIANGLE;
	this->number_boundaries = 3;

    this->boundary_type = boundary_type;
    this->neighbor_ID = neighbor_ID;

    this->nodal_coordinates = nodal_coordinates;

    this->number_bf = this->basis->GetNumberBasisFunctions();
    if (basis_geom != nullptr) {
        this->number_bf_geom = this->basis_geom->GetNumberBasisFunctions();
    }
    else {
        this->number_bf_geom = 3;
    }

    this->number_gp_internal = this->basis->GetIntegrationRuleInternal()->GetNumberGP();
    this->number_gp_boundary = this->basis->GetIntegrationRuleBoundary()->GetNumberGP();

	this->allocate_memory();

	for (int i = 0; i < this->number_boundaries; i++) {
		this->interfaces[i] = nullptr;
		this->interface_owner[i] = false;
	}

    //COMPUTE GEOMETRY AND M_INV
    if (this->basis_geom == nullptr) {
        this->orthogonal = this->basis->GetOrthogonal();
        this->m_inv = this->basis->GetMInv();

        double** J = new double*[2];
        J[0] = new double[2];
        J[1] = new double[2];

        J[0][0] = (this->nodal_coordinates[X][1] - this->nodal_coordinates[X][0]) / 2.0;
        J[0][1] = (this->nodal_coordinates[X][2] - this->nodal_coordinates[X][0]) / 2.0;
        J[1][0] = (this->nodal_coordinates[Y][1] - this->nodal_coordinates[Y][0]) / 2.0;
        J[1][1] = (this->nodal_coordinates[Y][2] - this->nodal_coordinates[Y][0]) / 2.0;

        this->det_J_internal[0] = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        this->J_inv_t_internal[0][0][0] = J[1][1] / (this->det_J_internal[0]);
        this->J_inv_t_internal[0][1][0] = -J[1][0] / (this->det_J_internal[0]);
        this->J_inv_t_internal[1][0][0] = -J[0][1] / (this->det_J_internal[0]);
        this->J_inv_t_internal[1][1][0] = J[0][0] / (this->det_J_internal[0]);

        //printf("%f %f\n%f %f\n\n", this->J_inv_t_area[0][0][0], this->J_inv_t_area[0][1][0],
        //	this->J_inv_t_area[1][0][0], this->J_inv_t_area[1][1][0]);

        this->surface_J_boundary[0][0] = sqrt(pow((J[0][0] - J[0][1]), 2.0) + pow((J[1][0] - J[1][1]), 2.0));
        this->surface_J_boundary[1][0] = sqrt(pow(J[0][1], 2.0) + pow(J[1][1], 2.0));
        this->surface_J_boundary[2][0] = sqrt(pow(J[0][0], 2.0) + pow(J[1][0], 2.0));

        double cw = this->det_J_internal[0] / abs(this->det_J_internal[0]); //CW or CCW

        this->normal_boundary[X][0][0] = (J[1][1] - J[1][0]) / this->surface_J_boundary[0][0] * cw;
        this->normal_boundary[X][1][0] = -J[1][1] / this->surface_J_boundary[1][0] * cw;
        this->normal_boundary[X][2][0] = J[1][0] / this->surface_J_boundary[2][0] *cw;

        this->normal_boundary[Y][0][0] = (J[0][0] - J[0][1]) / this->surface_J_boundary[0][0] * cw;
        this->normal_boundary[Y][1][0] = J[0][1] / this->surface_J_boundary[1][0] * cw;
        this->normal_boundary[Y][2][0] = -J[0][0] / this->surface_J_boundary[2][0] * cw;

        delete[] J[0];
        delete[] J[1];

        delete[] J;
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
	double** phi = this->basis->GetPhiPostProcessorCell();
	double* u_post = new double[N_DIV*N_DIV];

	for (int i = 0; i < N_DIV*N_DIV; i++) {
		u_post[i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < N_DIV*N_DIV; j++) {
			u_post[j] += this->u[u_flag][i] * phi[i][j];
		}
	}

	for (int i = 0; i < N_DIV*N_DIV; i++) {
		cell_data.push_back(u_post[i]);
	}

	delete[] u_post;
}

void ELEMENT::WritePointDataVTKTriangle(std::vector<double>& point_data, int u_flag) {
	double** phi = this->basis->GetPhiPostProcessorPoint();
	double* u_post = new double[(N_DIV + 1)*(N_DIV + 2) / 2];

	for (int i = 0; i < (N_DIV + 1)*(N_DIV + 2) / 2; i++) {
		u_post[i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < (N_DIV + 1)*(N_DIV + 2) / 2; j++) {
			u_post[j] += this->u[u_flag][i] * phi[i][j];
		}
	}

	for (int i = 0; i < (N_DIV + 1)*(N_DIV + 2) / 2; i++) {
		point_data.push_back(u_post[i]);
	}

	delete[] u_post;
}