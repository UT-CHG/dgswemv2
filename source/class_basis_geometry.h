#ifndef CLASS_BASIS_GEOMETRY_H
#define CLASS_BASIS_GEOMETRY_H

#include "class_integration.h"

class BASIS_GEOM_2D {
private:
    int p_geom;

    INTEGRATION_1D* integration_rule_line;
    INTEGRATION_2D* integration_rule_area;
    
    int number_bf_geom;

    double** dN_dz1_area;
    double** dN_dz2_area;
    double*** dN_dz1_edge;
    double*** dN_dz2_edge;

public:
    BASIS_GEOM_2D(int, INTEGRATION_1D*, INTEGRATION_2D*);
    ~BASIS_GEOM_2D();

    int GetPolynomial();
    int GetNumberBasisFunctions();
};

#endif