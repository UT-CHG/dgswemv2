#include "../../general_definitions.hpp"

#include "../../class_mesh.hpp"
#include "../../stepper.hpp"

uint main(uint argc, const char* argv[]){
	MESH* mesh = new MESH(1);

	mesh->RectangularDomainTest(90000.0, 45000.0, 100, 2, TRIANGLE);

	mesh->Solve();

	delete mesh;
}