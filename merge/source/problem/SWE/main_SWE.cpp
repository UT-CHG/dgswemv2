#include "../../general_definitions.hpp"

#include "swe_data.hpp"
#include "swe_kernels.hpp"
#include "swe_boundary_conditions.hpp"
#include "swe_definitions.hpp"

#include "../../geometry/mesh.hpp"
#include "../../stepper.hpp"
#include "../../initialize_mesh.hpp"
#include "../../run_simulation.hpp"

int main(int argc, const char* argv[]) {
	Geometry::MeshType<SWE::Data, SWE::Land, SWE::Tidal> mesh(2);
	Geometry::initialize_mesh(mesh);

	Stepper stepper(2, 2, 1.);

	run_simulation(344000.0, stepper, mesh);
}