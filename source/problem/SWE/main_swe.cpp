#include "../../general_definitions.hpp"

#include "../../mesh_metadata.hpp"
#include "../../stepper.hpp"

#include "../../initialize_mesh.hpp"
#include "../../initialize_data.hpp"
#include "../../run_simulation.hpp"

#include "../../ADCIRC_reader/adcirc_format.hpp"
#include "swe_definitions.hpp"
#include "swe_kernels.hpp"

int main(int argc, const char* argv[]) {
	AdcircFormat adcirc_file("sample_fort.14");
	MeshMetaData mesh_data(adcirc_file);

	auto mesh = initialize_mesh<SWE::Problem>(3, mesh_data);

	initialize_data(*mesh, adcirc_file);

	Stepper stepper(2, 2, 1.);

	auto t1 = std::chrono::high_resolution_clock::now();
	run_simulation<SWE::Problem>(86400.0, stepper, *mesh);
	auto t2 = std::chrono::high_resolution_clock::now();

	std::cout << "Time Elapsed (in us): "
		<< std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "\n";

	delete mesh;
}