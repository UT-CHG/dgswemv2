#include <iostream>
#include <fstream>
using namespace std;

#include "class_basis.h"
#include "basis_functions/basis_functions.h"

BASIS_TRI::BASIS_TRI(int p, AREA_INTEGRATION* area_rule, LINE_INTEGRATION* line_rule) {
	this->p = p;

	this->integration_rule_area = area_rule;
	this->integration_rule_line = line_rule;

	this->Dubiner();
}

BASIS_TRI::~BASIS_TRI(){
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
	delete[] phi_edge;
}

int BASIS_TRI::GetPolynomial() { return this->p; }

void BASIS_TRI::Dubiner() {
	int number_gp_area = this->integration_rule_area->GetNumberGP();
	int number_gp_edge = this->integration_rule_line->GetNumberGP();
	
	this->number_bf = (this->p + 1)*(this->p + 2) / 2;

	this->phi_area = new double*[this->number_bf];
	this->dphi_dz1_area = new double*[this->number_bf];
	this->dphi_dz2_area = new double*[this->number_bf];
	this->phi_edge = new double**[3];

	double* n1 = new double[number_gp_area];
	double* n2 = new double[number_gp_area];

	for (int i = 0; i < number_gp_area; i++) {
		n1[i] = 2 * (1 + this->integration_rule_area->GetZ1(i)) / (1 - this->integration_rule_area->GetZ2(i)) - 1;
		n2[i] = this->integration_rule_area->GetZ2(i);
	}

	int m = 0;
	for (int i = 0; i <= this->p; i++) {
		for (int j = 0; j <= this->p - i; j++) {
			this->phi_area[m] = new double[number_gp_area];
			this->dphi_dz1_area[m] = new double[number_gp_area];
			this->dphi_dz2_area[m] = new double[number_gp_area];

			dubiner_phi(i, j, number_gp_area, n1, n2, this->phi_area[m]);
			dubiner_dphi(i, j, number_gp_area, n1, n2, this->dphi_dz1_area[m], this->dphi_dz2_area[m]);

			m = m + 1;
		}
	}

	delete[] n1;
	delete[] n2;

	n1 = new double[number_gp_edge];
	n2 = new double[number_gp_edge];
	
	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			for (int j = 0; j < number_gp_edge; j++) {
				n1[j] = 1;
				n2[j] = this->integration_rule_line->GetZ(j);
			}
		}
		else if (i == 1) {
			for (int j = 0; j < number_gp_edge; j++) {
				n1[j] = -1;
				n2[j] = -this->integration_rule_line->GetZ(j);
			}
		}
		else if (i == 2) {
			for (int j = 0; j < number_gp_edge; j++) {
				n1[j] = this->integration_rule_line->GetZ(j);
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

	/*
	ofstream myfile;
	ofstream myfile_d;

	myfile.open("plot.txt");
	myfile_d.open("plot_d.txt");

	for (int i = 0; i < number_gp_area; i++) {
		myfile_d << this->integration_rule_area->GetZ1(i) << '\t';
		myfile_d << this->integration_rule_area->GetZ2(i) << '\t';
		myfile_d << this->dphi_dz1_area[20][i] << '\n';
	}

	for (int i = 0; i < number_gp_area; i++) {
		myfile << this->integration_rule_area->GetZ1(i) << '\t';
		myfile << this->integration_rule_area->GetZ2(i) << '\t';
		myfile << this->phi_area[20][i] << '\n';
	}

	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			for (int j = 0; j < number_gp_edge; j++) {
				myfile << -this->integration_rule_line->GetZ(j) << '\t';
				myfile << this->integration_rule_line->GetZ(j) << '\t';
				myfile << this->phi_edge[i][20][j] << '\n';
			}
		}
		else if (i == 1) {
			for (int j = 0; j < number_gp_edge; j++) {
				myfile << -1 << '\t';
				myfile << -this->integration_rule_line->GetZ(j) << '\t';
				myfile << this->phi_edge[i][20][j] << '\n';
			}
		}
		else if (i == 2) {
			for (int j = 0; j < number_gp_edge; j++) {
				myfile << this->integration_rule_line->GetZ(j) << '\t';
				myfile << -1 << '\t';
				myfile << this->phi_edge[i][20][j] << '\n';
			}
		}
	}

	myfile.close();
	*/

	double* w = new double[number_gp_area];
	for (int i = 0; i < number_gp_area; i++) {
		w[i] = this->integration_rule_area->GetWeight(i);
	}

	dubiner_test(this->p, number_gp_area, phi_area, w);

	delete[] n1;
	delete[] n2;
	delete[] w;
}