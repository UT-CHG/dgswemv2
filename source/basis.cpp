#include <iostream>
using namespace std;

#include "basis.h"
#include "basis_functions/basis_functions.h"

BASIS_TRI::BASIS_TRI(int p){
	this->p = p;
	
	AREA_INTEGRATION area_rule(2*p);
	LINE_INTEGRATION line_rule(2*p);

	this->integration_rule_area = &area_rule;
	this->integration_rule_line = &line_rule;

	this->Dubiner();
}

BASIS_TRI::~BASIS_TRI(){
	for (int i=0; i<this->number_bf; i++){
		delete [] phi_area[i]; 
	}
	delete [] phi_area;
}

int BASIS_TRI::GetPolynomial(){return this->p;}

void BASIS_TRI::Dubiner(){
	this->number_gp = this->integration_rule_area->GetNumberGP();
	this->number_bf = (this->p+1)*(this->p+2)/2;

	this->phi_area = new double*[this->number_bf];

	double* n1 = new double[this->number_gp];
	double* n2 = new double[this->number_gp];

	for (int i=0; i<this->number_gp; i++){		
		n1[i] = 2*(1 + this->integration_rule_area->GetZ1(i))/(1 - this->integration_rule_area->GetZ2(i)) - 1;
		n2[i] = this->integration_rule_area->GetZ2(i);
	}

	int m = 0;
	for (int i=0; i<=this->p; i++){
		for (int j=0; j<=this->p-i;j++){
			phi_area[m] = new double [this->number_gp];

			dubiner_phi(i,j,this->number_gp,n1,n2,phi_area[m]);

			m = m + 1;
		}
	}

	double* w = new double[this->number_gp];
	for (int i=0; i<this->number_gp; i++){
		w[i] = this->integration_rule_area->GetWeight(i);
	}

	dubiner_test(this->p, this->number_gp, phi_area, w);

	delete [] n1;
	delete [] n2;
	delete [] w;
}
