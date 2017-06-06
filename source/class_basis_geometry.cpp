#include <iostream>

#include "class_basis_geometry.h"

BASIS_GEOM::BASIS_GEOM(int p_geom, INTEGRATION* boundary_rule, INTEGRATION* internal_rule) {
    this->p_geom = p_geom;

    this->integration_rule_boundary = boundary_rule;
    this->integration_rule_internal = internal_rule;
}

BASIS_GEOM::~BASIS_GEOM(){ }