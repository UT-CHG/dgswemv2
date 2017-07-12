#ifndef RUN_SIMULATION_HPP
#define RUN_SIMULATION_HPP

template<typename ProblemType>
void run_simulation(double time_end, Stepper& stepper, typename ProblemType::mesh_type& mesh) {
	//we write these gross looking wrapper functions to append the stepper in a way that allows us to keep the
	//the nice std::for_each notation without having to define stepper within each element

	auto volume_kernel = [&stepper](auto& elt) {
		ProblemType::volume_kernel(stepper, elt);
	};

	auto source_kernel = [&stepper](auto& elt) {
		ProblemType::source_kernel(stepper, elt);
	};

	auto interface_kernel = [&stepper](auto& intface) {
		ProblemType::interface_kernel(stepper, intface);
	};

	auto boundary_kernel = [&stepper](auto& bound) {
		ProblemType::boundary_kernel(stepper, bound);
	};

	auto update_kernel = [&stepper](auto& elt) {
		ProblemType::update_kernel(stepper, elt);
	};

	auto swap_states_kernel = [&stepper](auto& elt) {
		ProblemType::swap_states_kernel(stepper, elt);
	};

	/*auto scru_sol = [&stepper](auto& elt) {
	  scrutinize_solution(stepper, elt);
	};*/

	uint nsteps = (uint)std::ceil(time_end / stepper.get_dt());
	uint n_stages = stepper.get_num_stages();

	auto resize_data_container = [n_stages](auto& elt) {
		elt.data.resize(n_stages);
	};

	mesh.CallForEachElement(resize_data_container);

	ProblemType::write_VTK_data_kernel(stepper, mesh);

	for (uint step = 1; step <= nsteps; ++step) {
		for (uint stage = 0; stage < stepper.get_num_stages(); ++stage) {

			mesh.CallForEachElement(volume_kernel);

			mesh.CallForEachElement(source_kernel);

			mesh.CallForEachInterface(interface_kernel);

			mesh.CallForEachBoundary(boundary_kernel);

			mesh.CallForEachElement(update_kernel);

			++stepper;
		}

		mesh.CallForEachElement(swap_states_kernel);

		if (step % 300 == 0) {
			std::cout << "Step: " << step << "\n";
			ProblemType::write_VTK_data_kernel(stepper, mesh);
		}
	}
}

#endif