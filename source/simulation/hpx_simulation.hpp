#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    const InputParameters input;

    typename ProblemType::mesh_type mesh;
    Stepper stepper;

  public:
    HPXSimulation()
        : input(),
          mesh(input.polynomial_order, input.mesh_data._mesh_name),
          stepper(input.rk.nstages, input.rk.order, input.dt) {}
    HPXSimulation(std::string input_string, uint locality, uint thread)
        : input(input_string, locality, thread),
          mesh(input.polynomial_order, input.mesh_data._mesh_name),
          stepper(input.rk.nstages, input.rk.order, input.dt) {
        printf("Starting program with p=%d for %s mesh\n\n", input.polynomial_order, input.mesh_file_name.c_str());

        initialize_mesh<ProblemType>(this->mesh, input.mesh_data);
    }

    hpx::future<void> Run(double);
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);
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

    //ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    //ProblemType::write_modal_data_kernel(this->stepper, this->mesh);

    hpx::future<void> future = hpx::make_ready_future();

    for (uint step = 1; step <= nsteps; ++step) {
        for (uint stage = 0; stage < this->stepper.get_num_stages(); ++stage) {
     //       future =
      //          future.then([this, &volume_kernel, &source_kernel, &interface_kernel](hpx::future<void>&&) {
                                this->mesh.CallForEachElement(volume_kernel);

                                this->mesh.CallForEachElement(source_kernel);

                                this->mesh.CallForEachInterface(interface_kernel);
        //                    })
          //          .then([this, &boundary_kernel, &update_kernel, &scrutinize_solution_kernel](hpx::future<void>&&) {
                         this->mesh.CallForEachBoundary(boundary_kernel);

                         this->mesh.CallForEachElement(update_kernel);

                         this->mesh.CallForEachElement(scrutinize_solution_kernel);

                         ++(this->stepper);
            //         });
        }

//        future = future.then([this, step, &swap_states_kernel](hpx::future<void>&&) {
            this->mesh.CallForEachElement(swap_states_kernel);

            if (step % 360 == 0) {
                std::cout << "Step: " << step << "\n";
                //ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
                //ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
            }
  //      });
    }

    return future;
}

template <typename ProblemType>
class HPXSimulationClient : hpx::components::client_base<HPXSimulationClient<ProblemType>, HPXSimulation<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationClient<ProblemType>, HPXSimulation<ProblemType>>;

  public:
    HPXSimulationClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> Run(double time_end) {
        using ActionType = typename HPXSimulation<ProblemType>::RunAction;
        return hpx::async<ActionType>(this->get_id(), time_end);
    }
};

#endif
