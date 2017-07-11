#include "../../general_definitions.hpp"

#include "../../ADCIRC_reader/adcirc_format.hpp"
#include "../../ADCIRC_reader/mesh_metadata.hpp"

#include "../../stepper.hpp"
#include "../../initialize_mesh.hpp"
#include "../../run_simulation.hpp"

#include "swe_definitions.hpp"
#include "swe_kernels.hpp"

int main(int argc, const char* argv[]) {
	AdcircFormat adcirc_file("sample_fort.14");
	MeshMetaData mesh_data(adcirc_file);
 	
	Stepper stepper(2, 2, 1.);

	auto mesh = initialize_mesh<SWE::Problem>(2, mesh_data);
	run_simulation<SWE::Problem>(345600.0, stepper, *mesh);
	
	delete mesh;
}