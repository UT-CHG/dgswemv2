#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "hpx_communicator.hpp"

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    const InputParameters input;

    Stepper stepper;
    typename ProblemType::mesh_type mesh;
    std::vector<HPXCommunicator<ProblemType>> communicators;

    std::string log_file_name;

  public:
    HPXSimulation()
        : input(),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order, input.mesh_data._mesh_name) {}
    HPXSimulation(const std::string& input_string, uint locality, uint thread)
        : input(input_string, locality, thread),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order, input.mesh_data._mesh_name) {
        this->log_file_name = "output/" + input.mesh_data._mesh_name + "_log";

        std::ofstream log_file(this->log_file_name, std::ofstream::out);

        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << input.mesh_file_name
                 << " mesh" << std::endl << std::endl;

        initialize_mesh<ProblemType>(this->mesh, input.mesh_data);

        this->SetUpCommunication(locality, thread);
    }

    hpx::future<void> Run();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);

  private:
    void SetUpCommunication(uint, uint);

    void Send(uint, uint);
    std::vector<hpx::future<uint>> Receive(uint);
};

template <typename ProblemType>
hpx::future<void> HPXSimulation<ProblemType>::Run() {
    // we write these gross looking wrapper functions to append the stepper in a way that allows us to keep the
    // the nice std::for_each notation without having to define stepper within each element
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    log_file << "Launching Simulation!" << std::endl << std::endl;

    uint nsteps = (uint)std::ceil(this->input.T_end / this->stepper.get_dt());
    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    ProblemType::write_modal_data_kernel(this->stepper, this->mesh);

   hpx::future<void> timestep_future = hpx::make_ready_future();
   uint timestamp = 0;

    for (uint step = 1; step <= nsteps; step++) {
        for (uint stage = 0; stage < n_stages; stage++) {     
            timestep_future = timestep_future.then([this, timestamp](hpx::future<void>&& timestep_future) {
                std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
                log_file << "Current (time, stage): (" << this->stepper.get_t_at_curr_stage() << ','
                         << this->stepper.get_stage() << ')' << std::endl;
#endif
                std::vector<hpx::future<uint>> receive_futures = this->Receive(timestamp);

                this->Send(this->stepper.get_t_at_curr_stage(), timestamp);
#ifdef VERBOSE
                log_file << "Starting work before receive" << std::endl;
#endif
                auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

                auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

                auto interface_kernel = [this](auto& intface) {
                    ProblemType::interface_kernel(this->stepper, intface);
                };

                this->mesh.CallForEachElement(volume_kernel);

                this->mesh.CallForEachElement(source_kernel);

                this->mesh.CallForEachInterface(interface_kernel);
#ifdef VERBOSE
                log_file << "Finished work before receive" << std::endl
                         << "Starting to wait on receive with timestamp: " << timestamp << std::endl;
#endif
                return when_all(receive_futures).then([this](
                    hpx::future<std::vector<hpx::future<uint>>>&& ready_messages) {
                    std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
                    for (auto& message : ready_messages.get()) {
                        log_file << "Received message: " << message.get() << std::endl;
                    }

                    log_file << "Starting work after receive" << std::endl;
#endif
                    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

                    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

                    auto scrutinize_solution_kernel = [this](auto& elt) {
                        ProblemType::scrutinize_solution_kernel(this->stepper, elt);
                    };

                    this->mesh.CallForEachBoundary(boundary_kernel);

                    this->mesh.CallForEachElement(update_kernel);

                    this->mesh.CallForEachElement(scrutinize_solution_kernel);

                    ++(this->stepper);
#ifdef VERBOSE
                    log_file << "Finished work after receive" << std::endl << std::endl;
#endif
                });
            });

            timestamp++;
        }

        timestep_future = timestep_future.then([this, step](hpx::future<void>&& timestep_future) {
            auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

            this->mesh.CallForEachElement(swap_states_kernel);

            if (step % 360 == 0) {
                std::ofstream log_file(this->log_file_name, std::ofstream::app);

                log_file << "Step: " << step << std::endl;

                ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
                ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
            }
        });
    }
    
    log_file << "Finished scheduling tasks!" << std::endl << std::endl;

    return timestep_future;
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
void HPXSimulation<ProblemType>::Send(uint message, uint timestamp) {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    for (HPXCommunicator<ProblemType>& communicator : this->communicators) {
        communicator.Send(message, timestamp);
#ifdef VERBOSE
        log_file << this->mesh.GetMeshName() << " sent message with time stamp: " << timestamp << std::endl;
#endif
    }
}

template <typename ProblemType>
std::vector<hpx::future<uint>> HPXSimulation<ProblemType>::Receive(uint timestamp) {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    std::vector<hpx::future<uint>> receive_futures;
    receive_futures.reserve(this->communicators.size());

    for (HPXCommunicator<ProblemType>& communicator : this->communicators) {
        receive_futures.push_back(communicator.Receive(timestamp));
#ifdef VERBOSE
        log_file << this->mesh.GetMeshName() << " posted receive for message with time stamp: " << timestamp
                 << std::endl;
#endif
    }

    return receive_futures;
}

template <typename ProblemType>
class HPXSimulationClient : hpx::components::client_base<HPXSimulationClient<ProblemType>, HPXSimulation<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationClient<ProblemType>, HPXSimulation<ProblemType>>;

  public:
    HPXSimulationClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> Run() {
        using ActionType = typename HPXSimulation<ProblemType>::RunAction;
        return hpx::async<ActionType>(this->get_id());
    }
};

#endif
