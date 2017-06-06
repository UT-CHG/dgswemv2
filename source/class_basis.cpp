#include <iostream>

#include "class_basis.h"

BASIS::BASIS(int type, int p, INTEGRATION* boundary_rule, INTEGRATION* internal_rule) {
    this->p = p;

    this->integration_rule_boundary = boundary_rule;
    this->integration_rule_internal = internal_rule;

    switch (type) {
    case DUBINER_2D: this->Dubiner2D(); break;
    default:
        printf("\n");
        printf("BASIS CONSTRUCTOR - Fatal error!\n");
        printf("Undefined basis type = %d\n", type);
        exit(1);
    }
}

BASIS::~BASIS(){
    for (int i = 0; i < this->number_bf; i++) {
        delete[] this->phi_internal[i];
		
		delete[] this->phi_postprocessor_cell[i];
		delete[] this->phi_postprocessor_point[i];
	}

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->number_bf; j++) {
			delete[] this->dphi_dz_internal[i][j];
		}
		delete[] this->dphi_dz_internal[i];
	}

	for (int i = 0; i < this->number_boundaries; i++) {
		for (int j = 0; j < this->number_bf; j++) {
			delete[] this->phi_boundary[i][j];
		}
		delete[] this->phi_boundary[i];
	}

    delete[] this->phi_internal;
	delete[] this->dphi_dz_internal;
	delete[] this->phi_boundary;

	delete[] this->phi_postprocessor_cell;
	delete[] this->phi_postprocessor_point;

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

void BASIS::allocate_memory(int number_gp_internal, int number_gp_boundary) {
	this->phi_internal = new double*[this->number_bf];
	this->dphi_dz_internal = new double**[this->dimension];
	this->phi_boundary = new double**[this->number_boundaries];

	this->phi_postprocessor_cell = new double*[this->number_bf];
	this->phi_postprocessor_point = new double*[this->number_bf];

	for (int i = 0; i < this->number_bf; i++) {
		this->phi_internal[i] = new double[number_gp_internal];
		this->phi_postprocessor_cell[i] = new double[N_DIV*N_DIV];
		this->phi_postprocessor_point[i] = new double[(N_DIV + 1)*(N_DIV + 2) / 2];
	}

	for (int i = 0; i < this->dimension; i++) {
		this->dphi_dz_internal[i] = new double*[this->number_bf];
		for (int j = 0; j < this->number_bf; j++) {
			this->dphi_dz_internal[i][j] = new double[number_gp_internal];
		}
	}

	for (int i = 0; i < this->number_boundaries; i++) {
		this->phi_boundary[i] = new double*[this->number_bf];
		for (int j = 0; j < this->number_bf; j++) {
			this->phi_boundary[i][j] = new double[number_gp_boundary];
		}
	}

	if (this->orthogonal) {
		this->m_inv = new double*[1];
		this->m_inv[0] = new double[this->number_bf];
	}
	else if (!(this->orthogonal)) {
		this->m_inv = new double*[this->number_bf];
		for (int i = 0; i < this->number_bf; i++) {
			m_inv[i] = new double[this->number_bf];
		}
	}
}

void BASIS::Dubiner2D() {
    this->orthogonal = true;
	this->dimension = 2;
	this->number_boundaries = 3;

	this->number_bf = (this->p + 1)*(this->p + 2) / 2;

	int number_gp_internal = this->integration_rule_internal->GetNumberGP();
	int number_gp_boundary = this->integration_rule_boundary->GetNumberGP();

	this->allocate_memory(number_gp_internal, number_gp_boundary);

    double* n1 = new double[number_gp_internal];
    double* n2 = new double[number_gp_internal];

    double* z1 = this->integration_rule_internal->GetZ()[X];
    double* z2 = this->integration_rule_internal->GetZ()[Y];

    for (int i = 0; i < number_gp_internal; i++) {
        n1[i] = 2 * (1 + z1[i]) / (1 - z2[i]) - 1;
        n2[i] = z2[i];
    }

    int m = 0;
    for (int i = 0; i <= this->p; i++) {
        for (int j = 0; j <= this->p - i; j++) {
            dubiner_2d_phi(i, j, m, number_gp_internal, n1, n2, this->phi_internal);
            dubiner_2d_dphi(i, j, m, number_gp_internal, n1, n2, this->dphi_dz_internal);

            this->m_inv[0][m] = (2 * i + 1)*(i + j + 1) / 2.0;

            m = m + 1;
        }
    }

    delete[] n1;
    delete[] n2;

    n1 = new double[number_gp_boundary];
    n2 = new double[number_gp_boundary];
    
    double* z = this->integration_rule_boundary->GetZ()[X];

    for (int i = 0; i < this->number_boundaries; i++) {
        if (i == 0) {
            for (int j = 0; j < number_gp_boundary; j++) {
                n1[j] = 1;
                n2[j] = z[j];
            }
        }
        else if (i == 1) {
            for (int j = 0; j < number_gp_boundary; j++) {
                n1[j] = -1;
                n2[j] = -z[j];
            }
        }
        else if (i == 2) {
            for (int j = 0; j < number_gp_boundary; j++) {
                n1[j] = z[j];
                n2[j] = -1;
            }
        }

        m = 0;
        for (int j = 0; j <= this->p; j++) {
            for (int k = 0; k <= this->p - j; k++) {
                dubiner_2d_phi(j, k, m, number_gp_boundary, n1, n2, this->phi_boundary[i]);

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
            dubiner_2d_phi(i, j, m, N_DIV*N_DIV, n1, n2, this->phi_postprocessor_cell);

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

			dubiner_2d_phi(i, j, m, (N_DIV + 1)*(N_DIV + 2) / 2 - 1, n1, n2, this->phi_postprocessor_point);

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

	dubiner_2d_test(this->p, number_gp_internal, this->phi_internal, this->integration_rule_internal->GetWeight());
}