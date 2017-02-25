#include <iostream>
using namespace std;

#include "integration.h"
#include "integration_rules/integration_rules_area.h"
#include "integration_rules/integration_rules_line.h"

AREA_INTEGRATION::AREA_INTEGRATION(int p) {
	this->p = p;
	this->Dunavant();
}

AREA_INTEGRATION::~AREA_INTEGRATION() {
	delete[] w;
	delete[] z1;
	delete[] z2;
}

int AREA_INTEGRATION::GetPolynomial() { return this->p; }

int AREA_INTEGRATION::GetNumberGP() { return this->number_gp; }

double AREA_INTEGRATION::GetWeight(int n) { return this->w[n]; }

double AREA_INTEGRATION::GetZ1(int n) { return this->z1[n]; }

double AREA_INTEGRATION::GetZ2(int n) { return this->z2[n]; }

void AREA_INTEGRATION::Dunavant() {
	int polynomial = dunavant_degree(this->p);

	this->number_gp = dunavant_number_gp(polynomial);
	double* x = new double[number_gp];
	double* y = new double[number_gp];
	double* weight = new double[number_gp];

	dunavant_rule(polynomial, x, y, weight);

	/*for (int i = 0; i < number_gp; i++){
		cout << weight[i] << '\t';
		cout << x[i] << '\t';
		cout << y[i] << '\n';
	}*/

	this->w = new double[number_gp];
	this->z1 = new double[number_gp];
	this->z2 = new double[number_gp];

	for (int i = 0; i < number_gp; i++) {
		w[i] = 2 * weight[i];
		z1[i] = 2 * x[i] - 1;
		z2[i] = 2 * y[i] - 1;
	}

	dunavant_rule_test(polynomial, number_gp, z1, z2, w);

	/*for (int i = 0; i < number_gp; i++){
		cout << w[i] << '\t';
		cout << z1[i] << '\t';
		cout << z2[i] << '\n';
	}*/

	delete[] x;
	delete[] y;
	delete[] weight;
}

LINE_INTEGRATION::LINE_INTEGRATION(int p) {
	this->p = p;
	this->GaussLegendre();
}

LINE_INTEGRATION::~LINE_INTEGRATION() {
	delete[] w;
	delete[] z;
}

int LINE_INTEGRATION::GetPolynomial(){return this->p;}

int LINE_INTEGRATION::GetNumberGP(){return this->number_gp;}

double LINE_INTEGRATION::GetWeight(int n){return this->w[n];}

double LINE_INTEGRATION::GetZ(int n){return this->z[n];}

void LINE_INTEGRATION::GaussLegendre(){
	int polynomial = gausslegendre_degree(this->p);
	
	this->number_gp = gausslegendre_number_gp(polynomial);
	this->z = new double [number_gp];
	this->w = new double [number_gp];

	gausslegendre_rule(number_gp, z, w);

	gausslegendre_rule_test(polynomial, number_gp, z, w);

	/*for (int i = 0; i < number_gp; i++){
		cout << w[i] << '\t'; 
		cout << z[i] << '\n';
	}*/
}
