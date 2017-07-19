#include "../../general_definitions.hpp"

#include "../../mesh_metadata.hpp"
#include "../../stepper.hpp"

#include "../../initialize_mesh.hpp"
#include "../../initialize_data.hpp"
#include "../../run_simulation.hpp"

#include "../../ADCIRC_reader/adcirc_format.hpp"
#include "swe_problem.hpp"
#include "swe_kernels.hpp"

int main(int argc, const char* argv[]) {
	printf("Starting program %s with p=%d for %s mesh\n\n", argv[0], std::stoi(argv[1]), argv[2]);

	AdcircFormat adcirc_file(argv[2]);
	MeshMetaData mesh_data(adcirc_file);

	auto mesh = initialize_mesh<SWE::Problem>(std::stoi(argv[1]), mesh_data);

	initialize_data(*mesh, adcirc_file);

	Stepper stepper(2, 2, 1.);

	auto t1 = std::chrono::high_resolution_clock::now();
	run_simulation<SWE::Problem>(5*86400.0, stepper, *mesh);
	auto t2 = std::chrono::high_resolution_clock::now();

	std::cout << "Time Elapsed (in us): "
		<< std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "\n";

	delete mesh;
}