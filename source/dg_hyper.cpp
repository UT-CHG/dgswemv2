#include <iostream>
#include <fstream>
#include <iomanip>

#include "elements/class_element_tri.h"
#include "class_basis.h"
#include "class_basis_geometry.h"
//#include "class_mesh.h"

int main(int argc, const char* argv[])
{
    printf("\nHello World\n\n");

    int ID = 0;
    int p = 2;
    int p_geom = 1;

    double x[3] = { -1, 0, -1 };
    double y[3] = { -1, 0, 1 };

    INTEGRATION_2D area_rule(2 * p);
    INTEGRATION_1D line_rule(2 * p);
    BASIS_2D basis(p, &line_rule, &area_rule);
    BASIS_GEOM_2D basis_geom(p_geom, &line_rule, &area_rule);

    ELEMENT_TRI* element = new ELEMENT_TRI(ID, x, y, &basis);
    element->CreateInterfaces();

    element->~ELEMENT_TRI();
}
