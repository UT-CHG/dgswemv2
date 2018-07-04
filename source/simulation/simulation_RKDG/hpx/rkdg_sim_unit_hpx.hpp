#ifndef RKDG_SIM_UNIT_HPX_HPP
#define RKDG_SIM_UNIT_HPX_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "communication/hpx_communicator.hpp"

#include "simulation/writer.hpp"
#include "simulation/simulation_RKDG/load_balancer/base_model.hpp"
#include "simulation/simulation_RKDG/load_balancer/abstract_load_balancer_factory.hpp"

namespace RKDG {
template <typename ProblemType>
class HPXSimulationUnit
    : public hpx::components::migration_support<hpx::components::component_base<HPXSimulationUnit<ProblemType>>> {
  private:
    typename ProblemType::ProblemMeshType mesh;

    HPXCommunicator communicator;
    RKStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

    typename ProblemType::ProblemInputType problem_input;
    std::unique_ptr<LoadBalancer::SubmeshModel> submesh_model = nullptr;

  public:
    HPXSimulationUnit() = default;
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);
    HPXSimulationUnit(HPXSimulationUnit&& rhs) = default;

    HPXSimulationUnit& operator=(HPXSimulationUnit&& rhs) = default;

    hpx::future<void> FinishPreprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, FinishPreprocessor, FinishPreprocessorAction);

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Launch, LaunchAction);

    hpx::future<void> Step();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Step, StepAction);

    double ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, ResidualL2, ResidualL2Action);

    template <typename Archive>
    void save(Archive& ar, unsigned) const;

    template <typename Archive>
    void load(Archive& ar, unsigned);
    HPX_SERIALIZATION_SPLIT_MEMBER();

    // Do not rename this is overload member of the base class
    void on_migrated();

  private:
    hpx::future<void> Stage();

    hpx::future<void> Postprocessor();

    void SwapStates();
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

    this->problem_input = input.problem_input;
    this->submesh_model = LoadBalancer::AbstractFactory::create_submesh_model<ProblemType>(locality_id, submesh_id);

    assert(this->submesh_model);

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    ProblemType::initialize_problem_parameters(this->problem_input);

    ProblemType::preprocess_mesh_data(input);

    initialize_mesh<ProblemType, HPXCommunicator>(this->mesh, input, this->communicator, this->writer);

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

    this->mesh.CallForEachElement([n_stages](auto& elt) { elt.data.resize(n_stages + 1); });
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Stage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.GetTimeAtCurrentStage() << ','
                                  << this->stepper.GetStage() << ')' << std::endl;

        this->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    hpx::future<void> receive_future = this->communicator.ReceiveAll(this->stepper.GetTimestamp());

    this->mesh.CallForEachDistributedBoundary(
        [this](auto& dbound) { ProblemType::distributed_boundary_send_kernel(this->stepper, dbound); });

    this->communicator.SendAll(this->stepper.GetTimestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    if (this->parser.ParsingInput()) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Parsing input" << std::endl;
        }

        this->parser.ParseInput(this->stepper, this->mesh);
    }

    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); });

    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); });

    this->mesh.CallForEachInterface([this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); });

    this->mesh.CallForEachBoundary([this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); });

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work before receive" << std::endl
                                  << "Starting to wait on receive with timestamp: " << this->stepper.GetTimestamp()
                                  << std::endl;
    }

    return receive_future.then([this](auto&& f) {
        f.get();  // check for exceptions

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        this->mesh.CallForEachDistributedBoundary(
            [this](auto& dbound) { ProblemType::distributed_boundary_kernel(this->stepper, dbound); });

        this->mesh.CallForEachElement([this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); });

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

        ProblemType::postprocessor_parallel_post_receive_kernel(this->stepper, this->mesh);

        this->mesh.CallForEachElement([this](auto& elt) {
            bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

            if (nan_found)
                hpx::terminate();
        });

        ++(this->stepper);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished postprocessor work after receive" << std::endl << std::endl;
        }
    });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::SwapStates() {
    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); });

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Step() {
    hpx::future<void> step_future = hpx::make_ready_future();

    for (uint stage = 0; stage < this->stepper.GetNumStages(); stage++) {
        step_future = step_future.then([this](auto&& f) {
            f.get();
            return this->Stage();
        });

        step_future = step_future.then([this](auto&& f) {
            f.get();
            return this->Postprocessor();
        });
    }

    return step_future.then([this](auto&& f) {
        f.get();
        this->submesh_model->InStep(0, 0);
        this->SwapStates();
    });
}

template <typename ProblemType>
double HPXSimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    this->mesh.CallForEachElement([this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    });

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}

template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::save(Archive& ar, unsigned) const {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Departing from locality " << hpx::get_locality_id() << std::endl;
    }

    ar& stepper& writer& mesh& problem_input& parser& communicator& submesh_model;
}

template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::load(Archive& ar, unsigned) {
    ar& stepper& writer& mesh& problem_input& parser& communicator& submesh_model;

    this->writer.StartLog();

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Arriving on locality " << hpx::get_locality_id() << std::endl;
    }
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::on_migrated() {
    this->mesh.SetMasters();

    this->mesh.CallForEachElement([this](auto& elt) {
        using MasterType = typename std::remove_reference<decltype(elt)>::type::ElementMasterType;

        using MasterElementTypes = typename decltype(this->mesh)::MasterElementTypes;

        MasterType& master_elt =
            std::get<Utilities::index<MasterType, MasterElementTypes>::value>(this->mesh.GetMasters());

        elt.SetMaster(master_elt);

        elt.Initialize();
    });

    initialize_mesh_interfaces_boundaries<ProblemType, HPXCommunicator>(mesh, problem_input, communicator, writer);
}

template <typename ProblemType>
class HPXSimulationUnitClient
    : public hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>>;

  public:
    HPXSimulationUnitClient() = default;
    HPXSimulationUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}
    HPXSimulationUnitClient(hpx::id_type&& id) : BaseType(std::move(id)) {}

    static constexpr const char* GetBasename() { return "Simulation_Unit_Client_"; }

    hpx::future<void> FinishPreprocessor() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::FinishPreprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Launch() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::LaunchAction;
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
}
#endif