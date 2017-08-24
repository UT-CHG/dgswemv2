#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "communication/hpx_communicator.hpp"

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    InputParameters input;

    Stepper stepper;
    typename ProblemType::mesh_type mesh;
    HPXCommunicator comm;

    std::string log_file_name;

  public:
    HPXSimulation()
        : input(),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order) {}
    HPXSimulation(std::string input_string, uint locality, uint sbmsh_id)
        : input(input_string, locality, sbmsh_id),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order),
        comm(locality, sbmsh_id, std::string(input.mesh_file_name + "_meta")) {

        hpx::cout << "Hello from HPXSimulation constructor " << locality << " " << sbmsh_id << '\n';

        input.ReadMesh();
        {
            std::string& mesh_name = mesh.GetMeshName();
            mesh_name = input.mesh_data._mesh_name;
        }

        this->log_file_name = "output/" + input.mesh_data._mesh_name + "_log";

        std::ofstream log_file(this->log_file_name, std::ofstream::out);
        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << input.mesh_file_name
                 << " mesh\n\n";

        initialize_mesh<ProblemType, HPXCommunicator>(this->mesh, input.mesh_data, comm);
    }

    hpx::future<void> Run();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);

  private:
    void Send(uint, uint);
    std::vector<hpx::future<uint>> Receive(uint);
};

template <typename ProblemType>
hpx::future<void> HPXSimulation<ProblemType>::Run() {
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

    uint nsteps = (uint)std::ceil(this->input.T_end / this->stepper.get_dt());
    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    // ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    // ProblemType::write_modal_data_kernel(this->stepper, this->mesh);

    hpx::future<void> future = hpx::make_ready_future();

    uint timestamp = 0;
    for (uint step = 1; step <= nsteps; step++) {
        for (uint stage = 0; stage < this->stepper.get_num_stages(); stage++) {/*
            std::vector<hpx::future<uint>> receive_futures;

            future = future.then([this, timestamp, &receive_futures](auto&&) {
                this->Send(this->stepper.get_t_at_curr_stage(), timestamp);

                receive_futures = this->Receive(timestamp);
            });

            future = future.then([this, &volume_kernel, &source_kernel, &interface_kernel](auto&&) {
                std::ofstream log_file(this->log_file_name, std::ofstream::app);

                log_file << this->mesh.GetMeshName() << " starting work on kernels before receive\n";

                this->mesh.CallForEachElement(volume_kernel);

                this->mesh.CallForEachElement(source_kernel);

                this->mesh.CallForEachInterface(interface_kernel);

                log_file << this->mesh.GetMeshName() << " finished work on kernels before receive\n";
            });

            future = future.then([this, timestamp, &receive_futures](auto&&) {
                std::ofstream log_file(this->log_file_name, std::ofstream::app);

                log_file << this->mesh.GetMeshName() << " starting to wait on receive with timestamp: " << timestamp
                         << '\n';

                /* when_all(receive_futures)
                    .then([this](auto&& ready_messages) {
                         std::ofstream log_file(this->log_file_name, std::ofstream::app);

                         for (auto& message : ready_messages.get()) {
                             log_file << this->mesh.GetMeshName() << " received message: " << message.get() << '\n';
                         }
                     })
                     .get();
            });
                            for (uint& message : messages) {
                                log_file << this->mesh.GetMeshName() << " received message: " << message << '\n';
                            }
                        });
            
            future = future.then([this, &boundary_kernel, &update_kernel, &scrutinize_solution_kernel](auto&&) {
                std::ofstream log_file(this->log_file_name, std::ofstream::app);

                log_file << this->mesh.GetMeshName() << " starting work on kernels after receive\n";

                this->mesh.CallForEachBoundary(boundary_kernel);

                this->mesh.CallForEachElement(update_kernel);

                this->mesh.CallForEachElement(scrutinize_solution_kernel);

                ++(this->stepper);

                log_file << this->mesh.GetMeshName() << " finished work on kernels after receive\n";
            });

            timestamp++;*/
        }
/*
        future = future.then([this, step, &swap_states_kernel](hpx::future<void>&&) {
            this->mesh.CallForEachElement(swap_states_kernel);
            if (step % 360 == 0) {
                std::cout << "Step: " << step << "\n";

                // ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
                // ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
            }
        });*/
    }

    return future;
}

/*template <typename ProblemType>
void HPXSimulation<ProblemType>::Send(uint message, uint timestamp) {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    for (HPXCommunicator<ProblemType>& communicator : this->communicators) {
        communicator.Send(message, timestamp);
        log_file << this->mesh.GetMeshName() << " sent message with time stamp: " << timestamp << '\n';
    }
}

template <typename ProblemType>
std::vector<hpx::future<uint>> HPXSimulation<ProblemType>::Receive(uint timestamp) {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    std::vector<hpx::future<uint>> receive_futures;
    receive_futures.reserve(this->communicators.size());

    for (HPXCommunicator<ProblemType>& communicator : this->communicators) {
        receive_futures.push_back(communicator.Receive(timestamp));
        log_file << this->mesh.GetMeshName() << " posted receive for message with time stamp: " << timestamp << '\n';
    }

    return receive_futures;
    }*/

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
