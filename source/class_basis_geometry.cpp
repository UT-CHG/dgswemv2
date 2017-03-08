#include <iostream>

#include "class_basis_geometry.h"

BASIS_GEOM_2D::BASIS_GEOM_2D(int p_geom, INTEGRATION_1D* line_rule, INTEGRATION_2D* area_rule) {
    this->p_geom = p_geom;

    this->integration_rule_line = line_rule;
    this->integration_rule_area = area_rule;
}

BASIS_GEOM_2D::~BASIS_GEOM_2D(){ }

int BASIS_GEOM_2D::GetPolynomial() { return this->p_geom; }

int BASIS_GEOM_2D::GetNumberBasisFunctions() { return this->number_bf_geom; }