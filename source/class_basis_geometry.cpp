#include <iostream>
#include <fstream>

#include "class_basis_geometry.h"

BASIS_GEOM_TRI::BASIS_GEOM_TRI(int p_geom, AREA_INTEGRATION* area_rule, LINE_INTEGRATION* line_rule) {
    this->p_geom = p_geom;

    this->integration_rule_area = area_rule;
    this->integration_rule_line = line_rule;
}

BASIS_GEOM_TRI::~BASIS_GEOM_TRI(){
}