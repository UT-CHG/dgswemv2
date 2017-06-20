#include <iostream>
#include <fstream>
#include <iomanip>

#include "../SWE/swe_data.hpp"
#include "../../class_mesh.h"

int main(int argc, const char* argv[]){
	MESH* mesh = new MESH(1);

	mesh->RectangularDomainTest(90000.0, 45000.0, 50, 4, TRIANGLE);

	delete mesh;
}