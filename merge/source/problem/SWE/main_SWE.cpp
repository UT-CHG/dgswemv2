#include "../../general_definitions.hpp"

#include "../../stepper.hpp"
#include "../../initialize_mesh.hpp"
#include "../../run_simulation.hpp"

#include "swe_definitions.hpp"
#include "swe_kernels.hpp"

int main(int argc, const char* argv[]) {
	Stepper stepper(2, 2, 1.);

	auto mesh = initialize_mesh<SWE::Problem>(2);
	run_simulation<SWE::Problem>(344000.0, stepper, *mesh);
	
	delete mesh;
}