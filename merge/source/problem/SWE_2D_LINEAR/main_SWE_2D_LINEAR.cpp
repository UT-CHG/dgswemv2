#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../class_mesh.h"
#include "../../stepper.hpp"

int main(int argc, const char* argv[]){
	MESH* mesh = new MESH(1);

	mesh->RectangularDomainTest(90000.0, 45000.0, 50, 4, TRIANGLE);

	delete mesh;

	Stepper stepper(2, 2, 10.);

	/*initialize_data(mesh, mesh_file, stepper);

	std::cout << "Number of interior edges: " << mesh.get_num_interior_edges() << "\n";
	std::cout << "Number of boundary edges: " << mesh.get_num_boundary_edges() << "\n";

	Writer::ParaviewWriter writer(2);
	std::string fname = "test.vtk"; std::string var_name = "Bathymetry";
	writer.write_output(fname, mesh, var_name,
		[](double ze, double qx, double qy, double bath) {
		return bath;
	});

	//run simulation
	auto t1 = std::chrono::high_resolution_clock::now();
	run_simulation(86400, stepper, mesh);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Time Elapsed (in us): "
		<< std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "\n";*/
}