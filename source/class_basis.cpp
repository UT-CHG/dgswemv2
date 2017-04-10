#include <iostream>
#include <fstream>

#include "class_basis.h"

BASIS_2D::BASIS_2D(int type, int p, INTEGRATION_1D* line_rule, INTEGRATION_2D* area_rule) {
    this->p = p;

    this->integration_rule_line = line_rule;
    this->integration_rule_area = area_rule;

    switch (type) {
    case DUBINER: this->Dubiner(); break;
    default:
        printf("\n");
        printf("BASIS_2D - Fatal error!\n");
        printf("Undefined basis type = %d\n", type);
        exit(1);
    }
}

BASIS_2D::~BASIS_2D(){
    for (int i = 0; i < this->number_bf; i++) {
        delete[] this->phi_area[i];
        delete[] this->dphi_dz1_area[i];
        delete[] this->dphi_dz2_area[i];
        delete[] this->phi_postprocessor_cell[i];
		delete[] this->phi_postprocessor_point[i];
	}
    delete[] this->phi_area;
    delete[] this->dphi_dz1_area;
    delete[] this->dphi_dz2_area;
    delete[] this->phi_postprocessor_cell;
	delete[] this->phi_postprocessor_point;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < this->number_bf; j++) {
            delete[] this->phi_edge[i][j];
        }
        delete[] this->phi_edge[i];
    }
    delete[] this->phi_edge;

    if (this->orthogonal) {
        delete[] this->m_inv[0];
    }
    else if (!(this->orthogonal)) {
        for (int i = 0; i < this->number_bf; i++) {
            delete[] this->m_inv[i];
        }
    }
    delete[] this->m_inv;
}

void BASIS_2D::Dubiner() {
    this->orthogonal = true;

    int number_gp_area = this->integration_rule_area->GetNumberGP();
    int number_gp_edge = this->integration_rule_line->GetNumberGP();
    
    this->number_bf = (this->p + 1)*(this->p + 2) / 2;

    this->phi_area = new double*[this->number_bf];
    this->dphi_dz1_area = new double*[this->number_bf];
    this->dphi_dz2_area = new double*[this->number_bf];
    this->phi_edge = new double**[3];
    
    this->m_inv = new double*[1];
    this->m_inv[0] = new double[this->number_bf];

    this->phi_postprocessor_cell = new double*[this->number_bf];
	this->phi_postprocessor_point = new double*[this->number_bf];

    double* n1 = new double[number_gp_area];
    double* n2 = new double[number_gp_area];

    double* z1 = this->integration_rule_area->GetZ1();
    double* z2 = this->integration_rule_area->GetZ2();

    for (int i = 0; i < number_gp_area; i++) {
        n1[i] = 2 * (1 + z1[i]) / (1 - z2[i]) - 1;
        n2[i] = z2[i];
    }

    int m = 0;
    for (int i = 0; i <= this->p; i++) {
        for (int j = 0; j <= this->p - i; j++) {
            this->phi_area[m] = new double[number_gp_area];
            this->dphi_dz1_area[m] = new double[number_gp_area];
            this->dphi_dz2_area[m] = new double[number_gp_area];

            dubiner_phi(i, j, number_gp_area, n1, n2, this->phi_area[m]);
            dubiner_dphi(i, j, number_gp_area, n1, n2, this->dphi_dz1_area[m], this->dphi_dz2_area[m]);

            this->m_inv[0][m] = (2 * i + 1)*(i + j + 1) / 2.0;

            m = m + 1;
        }
    }

    delete[] n1;
    delete[] n2;

    n1 = new double[number_gp_edge];
    n2 = new double[number_gp_edge];
    
    double* z = this->integration_rule_line->GetZ();

    for (int i = 0; i < 3; i++) {
        if (i == 0) {
            for (int j = 0; j < number_gp_edge; j++) {
                n1[j] = 1;
                n2[j] = z[j];
            }
        }
        else if (i == 1) {
            for (int j = 0; j < number_gp_edge; j++) {
                n1[j] = -1;
                n2[j] = -z[j];
            }
        }
        else if (i == 2) {
            for (int j = 0; j < number_gp_edge; j++) {
                n1[j] = z[j];
                n2[j] = -1;
            }
        }

        this->phi_edge[i] = new double*[this->number_bf];

        m = 0;
        for (int j = 0; j <= this->p; j++) {
            for (int k = 0; k <= this->p - j; k++) {
                this->phi_edge[i][m] = new double[number_gp_edge];

                dubiner_phi(j, k, number_gp_edge, n1, n2, this->phi_edge[i][m]);

                m = m + 1;
            }
        }
    }

    delete[] n1;
    delete[] n2;

    n1 = new double[N_DIV*N_DIV];
    n2 = new double[N_DIV*N_DIV];

    z1 = new double[N_DIV*N_DIV];
    z2 = new double[N_DIV*N_DIV];

    double dz = 2.0 / N_DIV;

    int n_pt = 0;
    for (int i = 0; i < N_DIV; i++) {
        for (int j = 0; j < N_DIV - i; j++) {
            z1[n_pt] = -1.0 + dz*j + dz / 3.0; //CENTROID
            z2[n_pt] = -1.0 + dz*i + dz / 3.0;

            n1[n_pt] = 2 * (1 + z1[n_pt]) / (1 - z2[n_pt]) - 1;
            n2[n_pt] = z2[n_pt];

            n_pt = n_pt + 1;
        }
    }

    for (int i = 1; i < N_DIV; i++) {
        for (int j = 0; j < N_DIV - i; j++) {
            z1[n_pt] = -1.0 + dz*j + 2 * dz / 3.0; //CENTROID
            z2[n_pt] = -1.0 + dz*i - dz / 3.0; 

            n1[n_pt] = 2 * (1 + z1[n_pt]) / (1 - z2[n_pt]) - 1;
            n2[n_pt] = z2[n_pt];

            n_pt = n_pt + 1;
        }
    }

    m = 0;
    for (int i = 0; i <= this->p; i++) {
        for (int j = 0; j <= this->p - i; j++) {
            this->phi_postprocessor_cell[m] = new double[N_DIV*N_DIV];

            dubiner_phi(i, j, N_DIV*N_DIV, n1, n2, this->phi_postprocessor_cell[m]);

            m = m + 1;
        }
    }

    delete[] n1;
    delete[] n2;

	delete[] z1;
	delete[] z2;

	n1 = new double[(N_DIV + 1)*(N_DIV + 2) / 2];
	n2 = new double[(N_DIV + 1)*(N_DIV + 2) / 2];

	z1 = new double[(N_DIV + 1)*(N_DIV + 2) / 2];
	z2 = new double[(N_DIV + 1)*(N_DIV + 2) / 2];

	n_pt = 0;
	for (int i = 0; i < N_DIV; i++) {
		for (int j = 0; j <= N_DIV - i; j++) {
			z1[n_pt] = -1.0 + dz*j;
			z2[n_pt] = -1.0 + dz*i;

			n1[n_pt] = 2 * (1 + z1[n_pt]) / (1 - z2[n_pt]) - 1;
			n2[n_pt] = z2[n_pt];

			n_pt = n_pt + 1;
		}
	}

	m = 0;
	for (int i = 0; i <= this->p; i++) {
		for (int j = 0; j <= this->p - i; j++) {
			this->phi_postprocessor_point[m] = new double[(N_DIV + 1)*(N_DIV + 2) / 2];

			dubiner_phi(i, j, (N_DIV + 1)*(N_DIV + 2) / 2 - 1, n1, n2, this->phi_postprocessor_point[m]);

			m = m + 1;
		}
	}

	//SINGULAR POINT

	double phi_singular[11][11] = {
		{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	};

	m = 0;
	for (int i = 0; i <= this->p; i++) {
		for (int j = 0; j <= this->p - i; j++) {
			this->phi_postprocessor_point[m][(N_DIV + 1)*(N_DIV + 2) / 2 - 1] = phi_singular[i][j];

			m = m + 1;
		}
	}

	delete[] n1;
	delete[] n2;

    delete[] z1;
    delete[] z2;

	dubiner_test(this->p, number_gp_area, this->phi_area, this->integration_rule_area->GetWeight());
}