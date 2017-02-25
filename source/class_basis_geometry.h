#ifndef CLASS_BASIS_GEOMETRY_H
#define CLASS_BASIS_GEOMETRY_H

#include "class_integration.h"

class BASIS_GEOM_TRI {
private:
	int p_geom;

	AREA_INTEGRATION* integration_rule_area;
	LINE_INTEGRATION* integration_rule_line;

	int number_bf_geom;

	double** dN_dz1_area;
	double** dN_dz2_area;
	double*** dN_dz1_edge;
	double*** dN_dz2_edge;

public:
	BASIS_GEOM_TRI(int, AREA_INTEGRATION*, LINE_INTEGRATION*);
	~BASIS_GEOM_TRI();
};

#endif