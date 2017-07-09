#ifndef RUN_SIMULATION_HPP
#define RUN_SIMULATION_HPP

#include "general_definitions.hpp"

template<typename MeshType>
void run_simulation(double time_end, Stepper& rk_stepper, MeshType& mesh) {
	//we write these gross looking wrapper functions to append the rk_stepper in a way that allows us to keep the
	//the nice std::for_each notation without having to define rk_stepper within each element

	auto volume_kernel = [&rk_stepper](auto& elt) {
		SWE::volume_kernel(rk_stepper, elt);
	};

	auto source_kernel = [&rk_stepper](auto& elt) {
		SWE::source_kernel(rk_stepper, elt);
	};

	auto interface_kernel = [&rk_stepper](auto& intface) {
		SWE::interface_kernel(rk_stepper, intface);
	};

	auto boundary_kernel = [&rk_stepper](auto& bound) {
		SWE::boundary_kernel(rk_stepper, bound);
	};

	auto update_kernel = [&rk_stepper](auto& elt) {
		SWE::update_kernel(rk_stepper, elt);
	};

	auto swap_states_kernel = [&rk_stepper](auto& elt) {
		SWE::swap_states_kernel(rk_stepper, elt);
	};

	/*auto scru_sol = [&rk_stepper](auto& elt) {
	  scrutinize_solution(rk_stepper, elt);
	};*/

	uint nsteps = std::ceil(time_end / rk_stepper.get_dt());
	uint n_stages = rk_stepper.get_num_stages();

	auto resize_data_container = [n_stages](auto& elt) {
		elt.data.resize(n_stages + 1);
	};

	mesh.CallForEachElement(resize_data_container);

	SWE::write_VTK_data(rk_stepper, mesh);

	for (uint step = 1; step <= nsteps; ++step) {
		for (uint stage = 0; stage < rk_stepper.get_num_stages(); ++stage) {

			mesh.CallForEachElement(volume_kernel);

			//mesh.CallForEachElement(source_kernel);

			mesh.CallForEachInterface(interface_kernel);

			mesh.CallForEachBoundary(boundary_kernel);

			mesh.CallForEachElement(update_kernel);

			++rk_stepper;
		}

		mesh.CallForEachElement(swap_states_kernel);

		if (step % 300 == 0) {
			std::cout << "Step: " << step << "\n";
			SWE::write_VTK_data(rk_stepper, mesh);
		}
	}
}

#endif