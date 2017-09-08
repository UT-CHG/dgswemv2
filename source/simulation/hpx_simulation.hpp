#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "communication/hpx_communicator.hpp"

#include "utilities/file_exists.hpp"

template <typename ProblemType>
struct SimulationUnit : public hpx::components::simple_component_base<SimulationUnit<ProblemType>> {
    InputParameters input;

    Stepper stepper;
    typename ProblemType::ProblemMeshType mesh;
    HPXCommunicator communicator;

    std::string log_file_name;

    SimulationUnit() : input(), stepper(input.rk.nstages, input.rk.order, input.dt), mesh(input.polynomial_order) {}
    SimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id)
        : input(input_string, locality_id, submesh_id),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order),
          communicator(std::string(input.mesh_file_name + "_meta"), locality_id, submesh_id) {
        input.ReadMesh();

        mesh.SetMeshName(input.mesh_data.mesh_name);

        this->log_file_name = "output/" + input.mesh_data.mesh_name + "_log";

        std::ofstream log_file(this->log_file_name, std::ofstream::out);

        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << input.mesh_file_name
                 << " mesh" << std::endl << std::endl;

        initialize_mesh<ProblemType, HPXCommunicator>(this->mesh, input.mesh_data, communicator);
    }

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(SimulationUnit, Launch, LaunchAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(SimulationUnit, Stage, StageAction);

    void Timestep();
    HPX_DEFINE_COMPONENT_ACTION(SimulationUnit, Timestep, TimestepAction);
};

template <typename ProblemType>
void SimulationUnit<ProblemType>::Launch() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    log_file << std::endl << "Launching Simulation!" << std::endl << std::endl;

    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
}

template <typename ProblemType>
hpx::future<void> SimulationUnit<ProblemType>::Stage() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
    log_file << "Current (time, stage): (" << this->stepper.get_t_at_curr_stage() << ',' << this->stepper.get_stage()
             << ')' << std::endl;

    log_file << "Starting work before receive" << std::endl;
#endif
    auto distributed_boundary_send_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_send_kernel(this->stepper, dbound);
    };

    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    hpx::future<void> receive_future = this->communicator.ReceiveAll(this->stepper.get_timestamp());

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_send_kernel);

    this->communicator.SendAll(this->stepper.get_timestamp());

    this->mesh.CallForEachElement(volume_kernel);

    this->mesh.CallForEachElement(source_kernel);

    this->mesh.CallForEachInterface(interface_kernel);
#ifdef VERBOSE
    log_file << "Finished work before receive" << std::endl
             << "Starting to wait on receive with timestamp: " << this->stepper.get_timestamp() << std::endl;
#endif
    return receive_future.then([this](auto&&) {
        std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
        log_file << "Starting work after receive" << std::endl;
#endif
        auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

        auto distributed_boundary_kernel = [this](auto& dbound) {
            ProblemType::distributed_boundary_kernel(this->stepper, dbound);
        };

        auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

        auto scrutinize_solution_kernel = [this](auto& elt) {
            ProblemType::scrutinize_solution_kernel(this->stepper, elt);
        };

        this->mesh.CallForEachBoundary(boundary_kernel);

        this->mesh.CallForEachDistributedBoundary(distributed_boundary_kernel);

        this->mesh.CallForEachElement(update_kernel);

        this->mesh.CallForEachElement(scrutinize_solution_kernel);

        ++(this->stepper);
#ifdef VERBOSE
        log_file << "Finished work after receive" << std::endl << std::endl;
#endif
    });
}

template <typename ProblemType>
void SimulationUnit<ProblemType>::Timestep() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->stepper.get_step() % 360 == 0) {
        std::ofstream log_file(this->log_file_name, std::ofstream::app);

        log_file << "Step: " << this->stepper.get_step() << std::endl;

        ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
        ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
class SimulationUnitClient
    : hpx::components::client_base<SimulationUnitClient<ProblemType>, SimulationUnit<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<SimulationUnitClient<ProblemType>, SimulationUnit<ProblemType>>;

  public:
    SimulationUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> Launch() {
        using ActionType = typename SimulationUnit<ProblemType>::LaunchAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Stage() {
        using ActionType = typename SimulationUnit<ProblemType>::StageAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Timestep() {
        using ActionType = typename SimulationUnit<ProblemType>::TimestepAction;
        return hpx::async<ActionType>(this->get_id());
    }
};

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    uint n_steps;
    uint n_stages;

    std::vector<SimulationUnitClient<ProblemType>> simulation_unit_clients;

  public:
    HPXSimulation() = default;
    HPXSimulation(const std::string& input_string) {
        const uint locality_id = hpx::get_locality_id();
        const hpx::naming::id_type here = hpx::find_here();

        InputParameters input(input_string);

        this->n_steps = (uint)std::ceil(input.T_end / input.dt);
        this->n_stages = input.rk.nstages;

        std::string mesh_file_prefix = input.mesh_file_name;
        mesh_file_prefix.erase(mesh_file_prefix.size() - 3);
        mesh_file_prefix += '_' + std::to_string(locality_id) + '_';

        uint submesh_id = 0;
        while (Utilities::file_exists(mesh_file_prefix + std::to_string(submesh_id) + ".14")) {
            hpx::future<hpx::id_type> simulation_unit_id =
                hpx::new_<hpx::components::simple_component<SimulationUnit<ProblemType>>>(
                    here, input_string, locality_id, submesh_id);

            this->simulation_unit_clients.emplace_back(
                SimulationUnitClient<ProblemType>(std::move(simulation_unit_id)));

            ++submesh_id;
        }
    }

    hpx::future<void> Run();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);
};

template <typename ProblemType>
hpx::future<void> HPXSimulation<ProblemType>::Run() {
    std::vector<hpx::future<void>> simulation_futures;

    for (auto& sim_unit_client : this->simulation_unit_clients) {
        simulation_futures.push_back(sim_unit_client.Launch());
    }

    for (uint step = 1; step < this->n_steps; step++) {
        for (uint stage = 0; stage < this->n_stages; stage++) {
            for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
                simulation_futures[sim_id] =
                    simulation_futures[sim_id]
                        .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Stage(); });
            }
        }

        for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
            simulation_futures[sim_id] =
                simulation_futures[sim_id]
                    .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Timestep(); });
        }
    }

    return hpx::when_all(simulation_futures);
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
