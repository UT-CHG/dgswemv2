#ifndef CLASS_INTEGRATION_H
#define CLASS_INTEGRATION_H

#include "general_definitions.h"
#include "integration_rules/integration_rules_1D.h"
#include "integration_rules/integration_rules_2D.h"

class INTEGRATION {

private:
	int dimension;
    int p;
    
	int number_gp;
    double* w;
    double** z;

public:
    INTEGRATION(int, int);
    ~INTEGRATION();

	int GetPolynomial() { return this->p; }
	int GetNumberGP() { return this->number_gp; }
	double* GetWeight() { return this->w; }
	double** GetZ() { return this->z; }

private:
	void allocate_memory();

    void GaussLegendre1D();
	void Dunavant2D();
};
#endif