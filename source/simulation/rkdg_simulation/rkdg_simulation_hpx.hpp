#ifndef RKDG_SIMULATION_HPX_HPP
#define RKDG_SIMULATION_HPX_HPP

#include "../../general_definitions.hpp"

#include "../../preprocessor/input_parameters.hpp"
#include "../../preprocessor/initialize_mesh.hpp"
#include "../../communication/hpx_communicator.hpp"
#include "../../utilities/file_exists.hpp"

#include <hpx/util/unwrapped.hpp>

#include "../writer.hpp"

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
class RKDGSimulationHPXUnit : public hpx::components::simple_component_base<RKDGSimulationHPXUnit<ProblemType>> {
  private:
    typename ProblemType::ProblemMeshType mesh;

    HPXCommunicator communicator;
    RKDGStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

  public:
    RKDGSimulationHPXUnit() = default;
    RKDGSimulationHPXUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);

    hpx::future<void> FinishPreprocessor();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPXUnit, FinishPreprocessor, FinishPreprocessorAction);

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPXUnit, Launch, LaunchAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPXUnit, Stage, StageAction);

    hpx::future<void> Postprocessor();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPXUnit, Postprocessor, PostprocessorAction);

    void Step();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPXUnit, Step, StepAction);

    double ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPXUnit, ResidualL2, ResidualL2Action);
};

template <typename ProblemType>
RKDGSimulationHPXUnit<ProblemType>::RKDGSimulationHPXUnit(const std::string& input_string,
                                                  const uint locality_id,
                                                  const uint submesh_id) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string, locality_id, submesh_id);

    input.read_mesh();                         // read mesh meta data
    input.read_bcis();                         // read bc data
    input.read_dbmd(locality_id, submesh_id);  // read distributed boundary meta data

    this->mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);

    this->communicator = HPXCommunicator(input.mesh_input.dbmd_data);
    this->stepper      = RKDGStepper(input.stepper_input);
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

    ProblemType::initialize_data_parallel_pre_send_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);
}

template <typename ProblemType>
hpx::future<void> RKDGSimulationHPXUnit<ProblemType>::FinishPreprocessor() {
    hpx::future<void> receive_future = this->communicator.ReceivePreprocAll(this->stepper.GetTimestamp());

    this->communicator.SendPreprocAll(this->stepper.GetTimestamp());

    return receive_future.then(
        [this](auto&&) { ProblemType::initialize_data_parallel_post_receive_kernel(this->mesh); });
}

template <typename ProblemType>
void RKDGSimulationHPXUnit<ProblemType>::Launch() {
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
hpx::future<void> RKDGSimulationHPXUnit<ProblemType>::Stage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.GetTimeAtCurrentStage() << ','
                                  << this->stepper.GetStage() << ')' << std::endl
                                  << "Starting work before receive" << std::endl;
    }

    auto distributed_boundary_send_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_send_kernel(this->stepper, dbound);
    };

    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    hpx::future<void> receive_future = this->communicator.ReceiveAll(this->stepper.GetTimestamp());

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_send_kernel);

    this->communicator.SendAll(this->stepper.GetTimestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    if (this->parser.ParsingInput()) {
        this->parser.ParseInput(this->stepper, this->mesh);
    }

    this->mesh.CallForEachElement(volume_kernel);

    this->mesh.CallForEachElement(source_kernel);

    this->mesh.CallForEachInterface(interface_kernel);

    this->mesh.CallForEachBoundary(boundary_kernel);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work before receive" << std::endl
                                  << "Starting to wait on receive with timestamp: " << this->stepper.GetTimestamp()
                                  << std::endl;
    }

    return receive_future.then([this](auto&&) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        auto distributed_boundary_kernel = [this](auto& dbound) {
            ProblemType::distributed_boundary_kernel(this->stepper, dbound);
        };

        auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

        this->mesh.CallForEachDistributedBoundary(distributed_boundary_kernel);

        this->mesh.CallForEachElement(update_kernel);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    });
}

template <typename ProblemType>
hpx::future<void> RKDGSimulationHPXUnit<ProblemType>::Postprocessor() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging postprocessor data" << std::endl;
    }

    hpx::future<void> receive_future = this->communicator.ReceivePostprocAll(this->stepper.GetTimestamp());

    ProblemType::postprocessor_parallel_pre_send_kernel(this->stepper, this->mesh);

    this->communicator.SendPostprocAll(this->stepper.GetTimestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting postprocessor work before receive" << std::endl;
    }

    ProblemType::postprocessor_parallel_pre_receive_kernel(this->stepper, this->mesh);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished postprocessor work before receive" << std::endl
                                  << "Starting to wait on postprocessor receive with timestamp: "
                                  << this->stepper.GetTimestamp() << std::endl;
    }

    return receive_future.then([this](auto&&) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting postprocessor work after receive" << std::endl;
        }

        auto scrutinize_solution_kernel = [this](auto& elt) {
            bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

            if (nan_found)
                hpx::terminate();
        };

        ProblemType::postprocessor_parallel_post_receive_kernel(this->stepper, this->mesh);

        this->mesh.CallForEachElement(scrutinize_solution_kernel);

        ++(this->stepper);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished postprocessor work after receive" << std::endl << std::endl;
        }
    });
}

template <typename ProblemType>
void RKDGSimulationHPXUnit<ProblemType>::Step() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
double RKDGSimulationHPXUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachElement(compute_residual_L2_kernel);

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}

template <typename ProblemType>
class RKDGSimulationHPXUnitClient
    : hpx::components::client_base<RKDGSimulationHPXUnitClient<ProblemType>, RKDGSimulationHPXUnit<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<RKDGSimulationHPXUnitClient<ProblemType>, RKDGSimulationHPXUnit<ProblemType>>;

  public:
    RKDGSimulationHPXUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> FinishPreprocessor() {
        using ActionType = typename RKDGSimulationHPXUnit<ProblemType>::FinishPreprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Launch() {
        using ActionType = typename RKDGSimulationHPXUnit<ProblemType>::LaunchAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Stage() {
        using ActionType = typename RKDGSimulationHPXUnit<ProblemType>::StageAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Postprocessor() {
        using ActionType = typename RKDGSimulationHPXUnit<ProblemType>::PostprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Step() {
        using ActionType = typename RKDGSimulationHPXUnit<ProblemType>::StepAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<double> ResidualL2() {
        using ActionType = typename RKDGSimulationHPXUnit<ProblemType>::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }
};

template <typename ProblemType>
class RKDGSimulationHPX : public hpx::components::simple_component_base<RKDGSimulationHPX<ProblemType>> {
  private:
    uint n_steps;
    uint n_stages;

    std::vector<RKDGSimulationHPXUnitClient<ProblemType>> simulation_unit_clients;

  public:
    RKDGSimulationHPX() = default;
    RKDGSimulationHPX(const std::string& input_string);

    hpx::future<void> Run();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPX, Run, RunAction);

    hpx::future<double> ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(RKDGSimulationHPX, ResidualL2, ResidualL2Action);
};

template <typename ProblemType>
RKDGSimulationHPX<ProblemType>::RKDGSimulationHPX(const std::string& input_string) {
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
            hpx::new_<hpx::components::simple_component<RKDGSimulationHPXUnit<ProblemType>>>(
                here, input_string, locality_id, submesh_id);

        this->simulation_unit_clients.emplace_back(std::move(simulation_unit_id));

        ++submesh_id;
    }
}

template <typename ProblemType>
hpx::future<void> RKDGSimulationHPX<ProblemType>::Run() {
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

            for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
                simulation_futures[sim_id] = simulation_futures[sim_id].then(
                    [this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Postprocessor(); });
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
hpx::future<double> RKDGSimulationHPX<ProblemType>::ResidualL2() {
    return ComputeL2Residual(this->simulation_unit_clients);
}

template <typename ProblemType>
class RKDGSimulationHPXClient : hpx::components::client_base<RKDGSimulationHPXClient<ProblemType>, RKDGSimulationHPX<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<RKDGSimulationHPXClient<ProblemType>, RKDGSimulationHPX<ProblemType>>;

  public:
    RKDGSimulationHPXClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}

    hpx::future<void> Run() {
        using ActionType = typename RKDGSimulationHPX<ProblemType>::RunAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<double> ResidualL2() {
        using ActionType = typename RKDGSimulationHPX<ProblemType>::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }
};

#endif
