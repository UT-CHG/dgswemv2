#ifndef CLASS_BASIS_GEOMETRY_H
#define CLASS_BASIS_GEOMETRY_H

#include "class_integration.h"

class BASIS_GEOM {
private:
	int dimension;
	int number_boundaries;

    INTEGRATION* integration_rule_boundary;
    INTEGRATION* integration_rule_internal;
    
    int p_geom;
    int number_bf_geom;

    double*** dN_dz_internal;
    double**** dN_dz1_boundary;

public:
    BASIS_GEOM(int, INTEGRATION*, INTEGRATION*);
    ~BASIS_GEOM();

	int GetPolynomial() { return this->p_geom; }
	int GetNumberBasisFunctions() { return this->number_bf_geom; }
};

#endif