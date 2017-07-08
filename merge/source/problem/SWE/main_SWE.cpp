#include "../../general_definitions.hpp"

#include "swe_data.hpp"
#include "swe_kernels.hpp"
#include "swe_boundary_conditions.hpp"
#include "swe_definitions.hpp"

#include "../../geometry/mesh.hpp"
#include "../../stepper.hpp"
#include "../../initialize_mesh.hpp"
#include "../../run_simulation.hpp"

//#include "../../ADCIRC_reader/adcirc_format.hpp"

int main(int argc, const char* argv[]){
	/*MESH* mesh = new MESH(2);

	mesh->RectangularDomainTest(90000.0, 45000.0, 20, 10, TRIANGLE);

	mesh->Solve();*/

	Geometry::MeshType<SWE::Data, SWE::Land, SWE::Tidal> mesh(2);
	Geometry::initialize_mesh(mesh);

	Stepper stepper(2, 2, 1.);

	run_simulation(344000.0, stepper,mesh);


	//delete mesh;
}