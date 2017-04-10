#include <iostream>

#include "class_integration.h"

INTEGRATION_1D::INTEGRATION_1D(int type, int p) {
    this->p = p;

    switch (type) {
    case GAUSS_LEGENDRE: this->GaussLegendre(); break;
    default:
        printf("\n");
        printf("INTEGRATION_1D - Fatal error!\n");
        printf("Undefined line integraton type = %d\n", type);
        exit(1);
    }
}

INTEGRATION_1D::~INTEGRATION_1D() {
    delete[] this->w;
    delete[] this->z;
}

void INTEGRATION_1D::GaussLegendre() {
    int polynomial = gausslegendre_degree(this->p);

    this->number_gp = gausslegendre_number_gp(polynomial);
    this->z = new double[this->number_gp];
    this->w = new double[this->number_gp];

    gausslegendre_rule(this->number_gp, this->z, this->w);

    gausslegendre_rule_test(polynomial, this->number_gp, this->z, this->w);
}

INTEGRATION_2D::INTEGRATION_2D(int type, int p) {
    this->p = p;

    switch (type) {
    case DUNAVANT: this->Dunavant(); break;
    default:
        printf("\n");
        printf("INTEGRATION_2D - Fatal error!\n");
        printf("Undefined area integraton type = %d\n", type);
        exit(1);
    }
}

INTEGRATION_2D::~INTEGRATION_2D() {
    delete[] this->w;
    delete[] this->z1;
    delete[] this->z2;
}

void INTEGRATION_2D::Dunavant() {
    int polynomial = dunavant_degree(this->p);

    this->number_gp = dunavant_number_gp(polynomial);
    double* x = new double[this->number_gp];
    double* y = new double[this->number_gp];
    double* weight = new double[this->number_gp];

    dunavant_rule(polynomial, x, y, weight);

    this->w = new double[this->number_gp];
    this->z1 = new double[this->number_gp];
    this->z2 = new double[this->number_gp];

    for (int i = 0; i < this->number_gp; i++) {
        this->w[i] = 2 * weight[i];
        this->z1[i] = 2 * x[i] - 1;
        this->z2[i] = 2 * y[i] - 1;
    }

    dunavant_rule_test(polynomial, this->number_gp, this->z1, this->z2, this->w);

    delete[] x;
    delete[] y;
    delete[] weight;
}