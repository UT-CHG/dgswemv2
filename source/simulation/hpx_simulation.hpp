#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "hpx_communicator.hpp"

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    const InputParameters input;

    typename ProblemType::mesh_type mesh;
    Stepper stepper;

    std::vector<HPXCommunicator<ProblemType>> communicators;

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

        this->SetUpCommunication(locality, thread);
    }

    hpx::future<void> Run(double);
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);

  private:
    void SetUpCommunication(uint, uint);

    void Send(uint, double);
    std::vector<hpx::future<uint>> Receive(double);
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

    hpx::future<void> future = hpx::make_ready_future();
    std::vector<hpx::future<uint>> receive_futures;

    future = future.then([this, &resize_data_container](hpx::future<void>&&) {
        this->mesh.CallForEachElement(resize_data_container);

        // ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
        // ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
    });

    for (uint step = 1; step <= nsteps; ++step) {
        for (uint stage = 0; stage < this->stepper.get_num_stages(); ++stage) {
            future = future.then([this, &volume_kernel, &source_kernel, &interface_kernel](hpx::future<void>&&) {
                this->Send(this->stepper.get_stage(), this->stepper.get_t_at_curr_stage());

                this->mesh.CallForEachElement(volume_kernel);

                this->mesh.CallForEachElement(source_kernel);

                this->mesh.CallForEachInterface(interface_kernel);
            });

            when_all(this->Receive(this->stepper.get_t_at_curr_stage())).get();

            future =
                future.then([this, &boundary_kernel, &update_kernel, &scrutinize_solution_kernel](hpx::future<void>&&) {
                    this->mesh.CallForEachBoundary(boundary_kernel);

                    this->mesh.CallForEachElement(update_kernel);

                    this->mesh.CallForEachElement(scrutinize_solution_kernel);

                    ++(this->stepper);
                });
        }

        future = future.then([this, step, &swap_states_kernel](hpx::future<void>&&) {
            this->mesh.CallForEachElement(swap_states_kernel);

            if (step % 360 == 0) {
                std::cout << "Step: " << step << "\n";
                // ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
                // ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
            }
        });
    }

    return future;
}

template <typename ProblemType>
void HPXSimulation<ProblemType>::SetUpCommunication(uint locality, uint thread) {
    std::string in_location;
    std::vector<std::string> out_locations;

    in_location = std::to_string(locality) + '_' + std::to_string(thread);

    // Here we need to find all out localities based on mesh info

    // This is just to communicate between everyone(complete graph)
    const uint n_localities = hpx::find_all_localities().size();
    const uint n_threads = hpx::get_os_thread_count();

    out_locations.reserve(n_localities * n_threads);

    for (uint locality = 0; locality < n_localities; locality++) {
        for (uint thread = 0; thread < n_threads; thread++) {
            out_locations.push_back(std::to_string(locality) + '_' + std::to_string(thread));
        }
    }

    out_locations.erase(std::remove(out_locations.begin(), out_locations.end(), in_location), out_locations.end());

    this->communicators.reserve(out_locations.size());

    for (std::string& out_location : out_locations) {
        this->communicators.emplace_back(in_location, out_location);
    }
}

template <typename ProblemType>
void HPXSimulation<ProblemType>::Send(uint message, double timestamp) {
    for (HPXCommunicator<ProblemType>& communicator : this->communicators) {
        communicator.Send(message, timestamp);
    }    
}

template <typename ProblemType>
std::vector<hpx::future<uint>> HPXSimulation<ProblemType>::Receive(double timestamp) {
    std::vector<hpx::future<uint>> receive_futures(communicators.size());
    
    for (HPXCommunicator<ProblemType>& communicator : this->communicators) {
        receive_futures.push_back(communicator.Receive(timestamp));
    }

    return receive_futures;    
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
