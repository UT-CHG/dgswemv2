#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "communication/hpx_communicator.hpp"

#include "utilities/file_exists.hpp"

template <typename ProblemType>
class HPXSimulationUnit : public hpx::components::simple_component_base<HPXSimulationUnit<ProblemType>> {
  private:
    InputParameters input;

    Stepper stepper;
    typename ProblemType::ProblemMeshType mesh;
    HPXCommunicator communicator;

    std::string log_file_name;

  public:
    HPXSimulationUnit() : input(), stepper(input.rk.nstages, input.rk.order, input.dt), mesh(input.polynomial_order) {}
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id)
        : input(input_string, locality_id, submesh_id),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order),
          communicator(input.mesh_file_path.substr(0, input.mesh_file_path.find_last_of('.')) + ".dbmd",
                       locality_id,
                       submesh_id) {
        input.ReadMesh();

        mesh.SetMeshName(input.mesh_data.mesh_name);

        this->log_file_name = "output/" + input.mesh_data.mesh_name + "_log";

        std::ofstream log_file(this->log_file_name, std::ofstream::out);

        if (!log_file) {
            std::cerr << "Error in opening log file, presumably the output directory does not exists.\n";
        }

        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << mesh.GetMeshName() << " mesh"
                 << std::endl << std::endl;

        initialize_mesh<ProblemType, HPXCommunicator>(this->mesh, input.mesh_data, communicator);
    }

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Launch, LaunchAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Stage, StageAction);

    void Step();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Step, StepAction);
};

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Launch() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    log_file << std::endl << "Launching Simulation!" << std::endl << std::endl;

    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Stage() {
#ifdef VERBOSE
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

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
#ifdef VERBOSE
        std::ofstream log_file(this->log_file_name, std::ofstream::app);

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
void HPXSimulationUnit<ProblemType>::Step() {
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
class HPXSimulationUnitClient
    : hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>>;

  public:
    HPXSimulationUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> Launch() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::LaunchAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Stage() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::StageAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Step() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::StepAction;
        return hpx::async<ActionType>(this->get_id());
    }
};

template <typename ProblemType>
class HPXSimulation : public hpx::components::simple_component_base<HPXSimulation<ProblemType>> {
  private:
    uint n_steps;
    uint n_stages;

    std::vector<HPXSimulationUnitClient<ProblemType>> simulation_unit_clients;

  public:
    HPXSimulation() = default;
    HPXSimulation(const std::string& input_string) {
        const uint locality_id = hpx::get_locality_id();
        const hpx::naming::id_type here = hpx::find_here();

        InputParameters input(input_string);

        this->n_steps = (uint)std::ceil(input.T_end / input.dt);
        this->n_stages = input.rk.nstages;

        std::string submesh_file_prefix = input.mesh_file_path.substr(0, input.mesh_file_path.find_last_of('.')) + "_" +
                                          std::to_string(locality_id) + '_';
        std::string submesh_file_postfix =
            input.mesh_file_path.substr(input.mesh_file_path.find_last_of('.'), input.mesh_file_path.size());

        uint submesh_id = 0;

        while (Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix)) {
            hpx::future<hpx::id_type> simulation_unit_id =
                hpx::new_<hpx::components::simple_component<HPXSimulationUnit<ProblemType>>>(
                    here, input_string, locality_id, submesh_id);

            this->simulation_unit_clients.emplace_back(std::move(simulation_unit_id));

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

    for (uint step = 1; step <= this->n_steps; step++) {
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
                    .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Step(); });
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
