#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../class_mesh.h"
#include "../../stepper.hpp"
//#include "../../initialize_data.hpp"
//#include "../../run_simulation.hpp"

int main(int argc, const char* argv[]){
	MESH* mesh = new MESH(3);

	mesh->RectangularDomainTest(90000.0, 45000.0, 50, 4, TRIANGLE);

	mesh->Solve();

	delete mesh;
}