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
    }
    delete[] this->phi_area;
    delete[] this->dphi_dz1_area;
    delete[] this->dphi_dz2_area;

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

int BASIS_2D::GetPolynomial() { return this->p; }

int BASIS_2D::GetNumberBasisFunctions() { return this->number_bf; }

INTEGRATION_1D* BASIS_2D::GetIntegrationRuleLine() { return this->integration_rule_line; }

INTEGRATION_2D* BASIS_2D::GetIntegrationRuleArea() { return this->integration_rule_area; }

double** BASIS_2D::GetPhiArea() { return this->phi_area; };

double** BASIS_2D::GetDPhiDZ1Area() { return this->dphi_dz1_area; };

double** BASIS_2D::GetDPhiDZ2Area() { return this->dphi_dz2_area; };

double*** BASIS_2D::GetPhiEdge() { return this->phi_edge; };

bool BASIS_2D::GetOrthogonal() { return this->orthogonal; };

double** BASIS_2D::GetMInv() { return this->m_inv; };

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



    dubiner_test(this->p, number_gp_area, this->phi_area, this->integration_rule_area->GetWeight());

    delete[] n1;
    delete[] n2;
}
