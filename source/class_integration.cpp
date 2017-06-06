#include <iostream>

#include "class_integration.h"

INTEGRATION::INTEGRATION(int type, int p) {
    this->p = p;

    switch (type) {
    case GAUSS_LEGENDRE_1D: this->GaussLegendre1D(); break;
	case DUNAVANT_2D: this->Dunavant2D(); break;
    default:
        printf("\n");
        printf("INTEGRATION CONSTRUCTOR - Fatal error!\n");
		printf("Undefined numerical integraton rule = %d\n", type);
        exit(1);
    }
}

INTEGRATION::~INTEGRATION() {
	for (int i = 0; i < this->dimension; i++)
		delete[] this->z[i];

	delete[] this->w;
	delete[] this->z;
}

void INTEGRATION::allocate_memory() {
	this->w = new double[this->number_gp];
	this->z = new double*[this->dimension];

	for (int i = 0; i < this->dimension; i++)
		this->z[i] = new double[this->number_gp];
}

void INTEGRATION::GaussLegendre1D() {
	this->dimension = 1;
    this->p = gausslegendre_1d_degree(this->p);

    this->number_gp = gausslegendre_1d_number_gp(this->p);

	this->allocate_memory();

    gausslegendre_1d_rule(this->number_gp, this->z, this->w);

    gausslegendre_1d_rule_test(this->p, this->number_gp, this->z, this->w);
}

void INTEGRATION::Dunavant2D() {
	this->dimension = 2;
	this->p = dunavant_2d_degree(this->p);

    this->number_gp = dunavant_2d_number_gp(this->p);

	this->allocate_memory();

    dunavant_2d_rule(this->p, this->z, this->w);

    dunavant_2d_rule_test(this->p, this->number_gp, this->z, this->w);
}