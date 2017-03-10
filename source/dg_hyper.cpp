#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "elements/class_element_tri.h"
#include "class_basis.h"
#include "class_basis_geometry.h"
//#include "class_mesh.h"

int main(int argc, const char* argv[])
{
    int ID = 0;
    int p = 2;
    int p_geom = 1;

    double x[3] = { -1, 0, -1 };
    double y[3] = { -1, 0, 1 };

	INTEGRATION_1D* line_rule = new INTEGRATION_1D(2 * p);
	INTEGRATION_2D* area_rule = new INTEGRATION_2D(2 * p);
	
	BASIS_2D* basis = new BASIS_2D(p, line_rule, area_rule);
    BASIS_GEOM_2D* basis_geom = new BASIS_GEOM_2D(p_geom, line_rule, area_rule);

    ELEMENT_TRI* element = new ELEMENT_TRI(ID, x, y, basis);
    element->CreateInterfaces();

	delete element;

	delete basis;
	delete basis_geom;

	delete line_rule;
	delete area_rule;

	_CrtDumpMemoryLeaks();
}
