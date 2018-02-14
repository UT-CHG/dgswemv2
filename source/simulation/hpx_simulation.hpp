#ifndef HPX_SIMULATION_HPP
#define HPX_SIMULATION_HPP

#include "../general_definitions.hpp"

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "../communication/hpx_communicator.hpp"
#include "../utilities/file_exists.hpp"

#include <hpx/util/unwrapped.hpp>

#include "writer.hpp"

template <typename ClientType>
hpx::future<double> ComputeL2Residual(std::vector<ClientType>& clients) {
    std::vector<hpx::future<double>> res_futures;

    for (uint id = 0; id < clients.size(); id++) {
        res_futures.push_back(clients[id].ResidualL2());
    }

    return hpx::when_all(res_futures).then([](auto && res_futures)->double {
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
    InputParameters<typename ProblemType::ProblemInputType> input;

    Stepper stepper;
    
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;
    typename ProblemType::ProblemMeshType mesh;
    HPXCommunicator communicator;

  public:
    HPXSimulationUnit()
        : input(),
          stepper(this->input.rk.nstages, this->input.rk.order, this->input.dt),
          mesh(this->input.polynomial_order) {}
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id)
        : input(input_string, locality_id, submesh_id),
          stepper(this->input.rk.nstages, this->input.rk.order, this->input.dt),
          writer(input, locality_id, submesh_id),
          parser(input, locality_id, submesh_id),
          mesh(this->input.polynomial_order),
          communicator(this->input.mesh_file_name.substr(0, this->input.mesh_file_name.find_last_of('.')) + ".dbmd",
                       locality_id,
                       submesh_id) {
        ProblemType::initialize_problem_parameters(this->input.problem_input);

        this->input.ReadMesh();

        this->mesh.SetMeshName(this->input.mesh_data.mesh_name);

        if (this->writer.WritingLog()) {
            this->writer.StartLog();

            this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                      << input.mesh_data.mesh_name << " mesh" << std::endl << std::endl;
        }

        initialize_mesh<ProblemType, HPXCommunicator>(
            this->mesh, this->input.mesh_data, this->communicator, this->input.problem_input, this->writer);
    }

    hpx::future<void> Preprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Preprocessor, PreprocessorAction);

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Launch, LaunchAction);

    void Parse();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Parse, ParseAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Stage, StageAction);

    hpx::future<void> Postprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Postprocessor, PostprocessorAction);

    void Step();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Step, StepAction);

    double ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, ResidualL2, ResidualL2Action);
};

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Preprocessor() {
    hpx::future<void> receive_future = this->communicator.ReceivePreprocAll(this->stepper.get_timestamp());

    ProblemType::initialize_data_parallel_pre_send_kernel(this->mesh, this->input.mesh_data, this->input.problem_input);

    this->communicator.SendPreprocAll(this->stepper.get_timestamp());

    return receive_future.then([this](auto&&) {
        ProblemType::initialize_data_parallel_post_receive_kernel(
            this->mesh, this->input.mesh_data, this->input.problem_input);
    });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Launch() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages + 1); };

    this->mesh.CallForEachElement(resize_data_container);
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Parse() {
    if (this->parser.ParsingInput()) {
        this->parser.ParseInput(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Stage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.get_t_at_curr_stage() << ','
                                  << this->stepper.get_stage() << ')' << std::endl << "Starting work before receive"
                                  << std::endl;
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

    hpx::future<void> receive_future = this->communicator.ReceiveAll(this->stepper.get_timestamp());

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_send_kernel);

    this->communicator.SendAll(this->stepper.get_timestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    this->mesh.CallForEachElement(volume_kernel);

    this->mesh.CallForEachElement(source_kernel);

    this->mesh.CallForEachInterface(interface_kernel);

    this->mesh.CallForEachBoundary(boundary_kernel);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work before receive" << std::endl
                                  << "Starting to wait on receive with timestamp: " << this->stepper.get_timestamp()
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

        auto scrutinize_solution_kernel = [this](auto& elt) {
            ProblemType::scrutinize_solution_kernel(this->stepper, elt);
        };

        this->mesh.CallForEachDistributedBoundary(distributed_boundary_kernel);

        this->mesh.CallForEachElement(update_kernel);

        this->mesh.CallForEachElement(scrutinize_solution_kernel);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    });
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Postprocessor() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging postprocessor data" << std::endl;
    }

    hpx::future<void> receive_future = this->communicator.ReceivePostprocAll(this->stepper.get_timestamp());

    ProblemType::postprocessor_parallel_pre_send_kernel(this->stepper, this->mesh);

    this->communicator.SendPostprocAll(this->stepper.get_timestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting postprocessor work before receive" << std::endl;
    }

    ProblemType::postprocessor_parallel_pre_receive_kernel(this->stepper, this->mesh);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished postprocessor work before receive" << std::endl
                                  << "Starting to wait on postprocessor receive with timestamp: "
                                  << this->stepper.get_timestamp() << std::endl;
    }

    return receive_future.then([this](auto&&) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting postprocessor work after receive" << std::endl;
        }

        ProblemType::postprocessor_parallel_post_receive_kernel(this->stepper, this->mesh);

        ++(this->stepper);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished postprocessor work after receive" << std::endl << std::endl;
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

    hpx::future<void> Preprocessor() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::PreprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Launch() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::LaunchAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Parse() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::ParseAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Stage() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::StageAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Postprocessor() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::PostprocessorAction;
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
    HPXSimulation(const std::string& input_string) {
        const uint locality_id = hpx::get_locality_id();
        const hpx::naming::id_type here = hpx::find_here();

        InputParameters<typename ProblemType::ProblemInputType> input(input_string);

        this->n_steps = (uint)std::ceil(input.T_end / input.dt);
        this->n_stages = input.rk.nstages;

        std::string submesh_file_prefix = input.mesh_file_name.substr(0, input.mesh_file_name.find_last_of('.')) + "_" +
                                          std::to_string(locality_id) + '_';
        std::string submesh_file_postfix =
            input.mesh_file_name.substr(input.mesh_file_name.find_last_of('.'), input.mesh_file_name.size());

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

    hpx::future<double> ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, ResidualL2, ResidualL2Action);
};

template <typename ProblemType>
hpx::future<void> HPXSimulation<ProblemType>::Run() {
    std::vector<hpx::future<void>> simulation_futures;

    for (auto& sim_unit_client : this->simulation_unit_clients) {
        simulation_futures.push_back(sim_unit_client.Preprocessor());
    }

    for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
        simulation_futures[sim_id] =
            simulation_futures[sim_id]
                .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Launch(); });
    }

    for (uint step = 1; step <= this->n_steps; step++) {
        for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
            simulation_futures[sim_id] =
                simulation_futures[sim_id]
                    .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Parse(); });
        }

        for (uint stage = 0; stage < this->n_stages; stage++) {
            for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
                simulation_futures[sim_id] =
                    simulation_futures[sim_id]
                        .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Stage(); });
            }

            for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
                simulation_futures[sim_id] =
                    simulation_futures[sim_id]
                        .then([this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Postprocessor(); });
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

#endif
