#include "../../general_definitions.hpp"

#include "../../class_mesh.hpp"
#include "../../stepper.hpp"
#include "../../ADCIRC_reader/adcirc_format.hpp"

int main(int argc, const char* argv[]){
 	//Geometry::AdcircFormat mesh_file("sample_fort.14");

	//mesh_file.write_to("a.14"); 

	MESH* mesh = new MESH(2);

	mesh->RectangularDomainTest(90000.0, 45000.0, 20, 10, TRIANGLE);

	mesh->Solve();

	auto a = [](){ return 1; };
  	
	auto volume_kernel = [](auto& elt) {
    	//SWE::volume_kernel(elt);
	};

	delete mesh;
}