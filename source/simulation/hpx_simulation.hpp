#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    const InputParameters input;

    typename ProblemType::mesh_type mesh;
    Stepper stepper;

  public:
    HPXSimulation(std::string input_string, uint locality, uint thread)
        : input(input_string.c_str(), locality, thread),
          mesh(input.polynomial_order, input.mesh_data._mesh_name),
          stepper(input.rk.nstages, input.rk.order, input.dt) {
        printf("Starting program with p=%d for %s mesh\n\n", input.polynomial_order, input.mesh_file_name.c_str());

        initialize_mesh<ProblemType>(this->mesh, input.mesh_data);
    }

    hpx::future<void> Run(double);
};

template <typename ProblemType>
hpx::future<void> HPXSimulation<ProblemType>::Run(double time_end) {
    // we write these gross looking wrapper functions to append the stepper in a way that allows us to keep the
    // the nice std::for_each notation without having to define stepper within each element
    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    auto scrutinize_solution_kernel = [this](auto& elt) {
        ProblemType::scrutinize_solution_kernel(this->stepper, elt);
    };

    uint nsteps = (uint)std::ceil(time_end / this->stepper.get_dt());
    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    ProblemType::write_modal_data_kernel(this->stepper, this->mesh);

    hpx::future<void> future = hpx::make_ready_future();

    for (uint step = 1; step <= nsteps; ++step) {
        for (uint stage = 0; stage < this->stepper.get_num_stages(); ++stage) {
            this->mesh.CallForEachElement(volume_kernel);

            this->mesh.CallForEachElement(source_kernel);

            this->mesh.CallForEachInterface(interface_kernel);

            future = future.then([this, &boundary_kernel](hpx::future<void>&&) {
                this->mesh.CallForEachBoundary(boundary_kernel);
            });

            this->mesh.CallForEachBoundary(boundary_kernel);

            this->mesh.CallForEachElement(update_kernel);

            this->mesh.CallForEachElement(scrutinize_solution_kernel);

            ++(this->stepper);
        }

        this->mesh.CallForEachElement(swap_states_kernel);

        if (step % 360 == 0) {
            std::cout << "Step: " << step << "\n";
            ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
            ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
        }
    }

    return future;
}

#endif
