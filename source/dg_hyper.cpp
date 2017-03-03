#include <iostream>
#include <fstream>
#include <iomanip>

#include "class_element_tri.h"
#include "class_basis.h"
#include "class_basis_geometry.h"
#include "class_mesh.h"

int main( int argc, const char* argv[] )
{
    printf( "\nHello World\n\n" );

    int ID = 0;
    int p = 5;
    int p_geom = 1;

    double x[3] = { -1, 0, -1 };
    double y[3] = { -1, 0, 1 };

    AREA_INTEGRATION area_rule(2 * p);
    LINE_INTEGRATION line_rule(2 * p);
    BASIS_TRI basis(p, &area_rule, &line_rule);
    BASIS_GEOM_TRI basis_geom(p_geom, &area_rule, &line_rule);

    ELEMENT_TRI element(ID, p, p_geom, x, y, &basis);
}
