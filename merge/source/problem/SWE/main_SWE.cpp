#include "../../general_definitions.hpp"

#include "../../class_mesh.hpp"

#include "../../stepper.hpp"
#include "../../mesh.hpp"
#include "../../initialize_mesh.hpp"

#include "../../ADCIRC_reader/adcirc_format.hpp"

int main(int argc, const char* argv[]){
	/*MESH* mesh = new MESH(2);

	mesh->RectangularDomainTest(90000.0, 45000.0, 20, 10, TRIANGLE);

	mesh->Solve();*/

	Geometry::MeshType<SWE::Data> mesh;
	Geometry::initialize_mesh(1, mesh);

	//delete mesh;
}