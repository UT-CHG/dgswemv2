#ifndef EHDG_SIMULATION_HPX_HPP
#define EHDG_SIMULATION_HPX_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "communication/hpx_communicator.hpp"
#include "utilities/file_exists.hpp"

#include <hpx/util/unwrapped.hpp>

#include "simulation/writer.hpp"

namespace EHDG {
template <typename ClientType>
hpx::future<double> ComputeL2Residual(std::vector<ClientType>& clients) {
    std::vector<hpx::future<double>> res_futures;

    for (uint id = 0; id < clients.size(); id++) {
        res_futures.push_back(clients[id].ResidualL2());
    }

    return hpx::when_all(res_futures).then([](auto&& res_futures) -> double {
        std::vector<double> res = hpx::util::unwrap(res_futures.get());
        double combined_res{0};
        for (double r : res) {
            combined_res += r;
        }
        return combined_res;
    });
}

template <typename ProblemType>
class HPXSimulationUnit : public hpx::components::simple_component_base<HPXSimulationUnit<ProblemType>> {
  private:
    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemMeshSkeletonType mesh_skeleton;

    HPXCommunicator communicator;
    RKStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

  public:
    HPXSimulationUnit() = default;
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);

    hpx::future<void> FinishPreprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, FinishPreprocessor, FinishPreprocessorAction);

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Launch, LaunchAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Stage, StageAction);

    void Step();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Step, StepAction);

    double ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, ResidualL2, ResidualL2Action);
};

template <typename ProblemType>
HPXSimulationUnit<ProblemType>::HPXSimulationUnit(const std::string& input_string,
                                                  const uint locality_id,
                                                  const uint submesh_id) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string, locality_id, submesh_id);

    input.read_mesh();                         // read mesh meta data
    input.read_bcis();                         // read bc data
    input.read_dbmd(locality_id, submesh_id);  // read distributed boundary meta data

    this->mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);

    this->communicator = HPXCommunicator(input.mesh_input.dbmd_data);
    this->stepper      = RKStepper(input.stepper_input);
    this->writer       = Writer<ProblemType>(input.writer_input, locality_id, submesh_id);
    this->parser       = typename ProblemType::ProblemParserType(input, locality_id, submesh_id);

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    ProblemType::initialize_problem_parameters(input.problem_input);

    ProblemType::preprocess_mesh_data(input);

    initialize_mesh<ProblemType, HPXCommunicator>(this->mesh, input, this->communicator, this->writer);
    initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, this->writer);

    ProblemType::initialize_data_parallel_pre_send_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::FinishPreprocessor() {
    hpx::future<void> receive_future = this->communicator.ReceivePreprocAll(this->stepper.GetTimestamp());

    this->communicator.SendPreprocAll(this->stepper.GetTimestamp());

    return receive_future.then(
        [this](auto&&) { ProblemType::initialize_data_parallel_post_receive_kernel(this->mesh); });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Launch() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    uint n_stages = this->stepper.GetNumStages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages + 1); };

    this->mesh.CallForEachElement(resize_data_container);
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Stage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.GetTimeAtCurrentStage() << ','
                                  << this->stepper.GetStage() << ')' << std::endl;
    }

    auto global_distributed_boundary_kernel = [this](auto& dbound) {
        ProblemType::global_distributed_boundary_kernel(this->stepper, dbound);
    };

    auto global_interface_kernel = [this](auto& intface) {
        ProblemType::global_interface_kernel(this->stepper, intface);
    };

    auto global_boundary_kernel = [this](auto& bound) { ProblemType::global_boundary_kernel(this->stepper, bound); };

    auto global_edge_interface_kernel = [this](auto& edge_int) {
        ProblemType::global_edge_interface_kernel(this->stepper, edge_int);
    };

    auto global_edge_boundary_kernel = [this](auto& edge_bound) {
        ProblemType::global_edge_boundary_kernel(this->stepper, edge_bound);
    };

    auto local_volume_kernel = [this](auto& elt) { ProblemType::local_volume_kernel(this->stepper, elt); };

    auto local_source_kernel = [this](auto& elt) { ProblemType::local_source_kernel(this->stepper, elt); };

    auto local_interface_kernel = [this](auto& intface) {
        ProblemType::local_interface_kernel(this->stepper, intface);
    };

    auto local_boundary_kernel = [this](auto& bound) { ProblemType::local_boundary_kernel(this->stepper, bound); };

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    hpx::future<void> receive_future = this->communicator.ReceiveAll(this->stepper.GetTimestamp());

    this->mesh.CallForEachDistributedBoundary(global_distributed_boundary_kernel);

    this->communicator.SendAll(this->stepper.GetTimestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    if (this->parser.ParsingInput()) {
        this->parser.ParseInput(this->stepper, this->mesh);
    }

    /* Global Pre Receive Step */
    this->mesh.CallForEachInterface(global_interface_kernel);

    this->mesh.CallForEachBoundary(global_boundary_kernel);

    this->mesh_skeleton.CallForEachEdgeInterface(global_edge_interface_kernel);

    this->mesh_skeleton.CallForEachEdgeBoundary(global_edge_boundary_kernel);
    /* Global Pre Receive Step */

    /* Local Pre Receive Step */
    this->mesh.CallForEachElement(local_volume_kernel);

    this->mesh.CallForEachElement(local_source_kernel);

    this->mesh.CallForEachInterface(local_interface_kernel);

    this->mesh.CallForEachBoundary(local_boundary_kernel);
    /* Local Pre Receive Step */

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work before receive" << std::endl
                                  << "Starting to wait on receive with timestamp: " << this->stepper.GetTimestamp()
                                  << std::endl;
    }

    return receive_future.then([this](auto&&) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        auto global_edge_distributed_kernel = [this](auto& edge_dbound) {
            ProblemType::global_edge_distributed_kernel(this->stepper, edge_dbound);
        };

        auto local_distributed_boundary_kernel = [this](auto& dbound) {
            ProblemType::local_distributed_boundary_kernel(this->stepper, dbound);
        };

        auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

        auto scrutinize_solution_kernel = [this](auto& elt) {
            bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

            if (nan_found)
                hpx::terminate();
        };

        /* Global Post Receive Step */
        this->mesh_skeleton.CallForEachEdgeDistributed(global_edge_distributed_kernel);
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        this->mesh.CallForEachDistributedBoundary(local_distributed_boundary_kernel);

        this->mesh.CallForEachElement(update_kernel);

        this->mesh.CallForEachElement(scrutinize_solution_kernel);
        /* Local Post Receive Step */

        ++(this->stepper);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Step() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
double HPXSimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachElement(compute_residual_L2_kernel);

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}

template <typename ProblemType>
class HPXSimulationUnitClient
    : hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>>;

  public:
    HPXSimulationUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> FinishPreprocessor() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::FinishPreprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

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

    hpx::future<double> ResidualL2() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::ResidualL2Action;
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
    HPXSimulation(const std::string& input_string);

    hpx::future<void> Run();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);

    hpx::future<double> ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, ResidualL2, ResidualL2Action);
};

template <typename ProblemType>
HPXSimulation<ProblemType>::HPXSimulation(const std::string& input_string) {
    const uint locality_id          = hpx::get_locality_id();
    const hpx::naming::id_type here = hpx::find_here();

    InputParameters<typename ProblemType::ProblemInputType> input(input_string);

    this->n_steps  = (uint)std::ceil(input.stepper_input.run_time / input.stepper_input.dt);
    this->n_stages = input.stepper_input.nstages;

    std::string submesh_file_prefix =
        input.mesh_input.mesh_file_name.substr(0, input.mesh_input.mesh_file_name.find_last_of('.')) + "_" +
        std::to_string(locality_id) + '_';
    std::string submesh_file_postfix = input.mesh_input.mesh_file_name.substr(
        input.mesh_input.mesh_file_name.find_last_of('.'), input.mesh_input.mesh_file_name.size());

    uint submesh_id = 0;

    while (Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix)) {
        hpx::future<hpx::id_type> simulation_unit_id =
            hpx::new_<hpx::components::simple_component<HPXSimulationUnit<ProblemType>>>(
                here, input_string, locality_id, submesh_id);

        this->simulation_unit_clients.emplace_back(std::move(simulation_unit_id));

        ++submesh_id;
    }
}

template <typename ProblemType>
hpx::future<void> HPXSimulation<ProblemType>::Run() {
    std::vector<hpx::future<void>> simulation_futures;

    for (auto& sim_unit_client : this->simulation_unit_clients) {
        simulation_futures.push_back(sim_unit_client.FinishPreprocessor());
    }

    for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
        simulation_futures[sim_id] = simulation_futures[sim_id].then(
            [this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Launch(); });
    }

    for (uint step = 1; step <= this->n_steps; step++) {
        for (uint stage = 0; stage < this->n_stages; stage++) {
            for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
                simulation_futures[sim_id] = simulation_futures[sim_id].then(
                    [this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Stage(); });
            }
        }

        for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
            simulation_futures[sim_id] = simulation_futures[sim_id].then(
                [this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Step(); });
        }
    }
    return hpx::when_all(simulation_futures);
}

template <typename ProblemType>
hpx::future<double> HPXSimulation<ProblemType>::ResidualL2() {
    return ComputeL2Residual(this->simulation_unit_clients);
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

    hpx::future<double> ResidualL2() {
        using ActionType = typename HPXSimulation<ProblemType>::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }
};
}

#endif
