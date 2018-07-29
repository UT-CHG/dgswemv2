#ifndef SIM_UNIT_HPX_HPP
#define SIM_UNIT_HPX_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "communication/hpx_communicator.hpp"

#include "simulation/writer.hpp"
#include "simulation/hpx/load_balancer/base_model.hpp"
#include "simulation/hpx/load_balancer/abstract_load_balancer_factory.hpp"

template <typename ProblemType>
struct HPXSimulationUnit
    : public hpx::components::migration_support<hpx::components::component_base<HPXSimulationUnit<ProblemType>>> {
    // *** //
    typename ProblemType::ProblemDiscretizationType discretization;

    HPXCommunicator communicator;
    RKStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

    typename ProblemType::ProblemInputType problem_input;
    std::unique_ptr<LoadBalancer::SubmeshModel> submesh_model = nullptr;

    HPXSimulationUnit() = default;
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);
    HPXSimulationUnit(HPXSimulationUnit&& rhs) = default;

    HPXSimulationUnit& operator=(HPXSimulationUnit&& rhs) = default;

    hpx::future<void> FinishPreprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, FinishPreprocessor, FinishPreprocessorAction);

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Launch, LaunchAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Stage, StageAction);

    void SwapStates();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, SwapStates, SwapStatesAction);

    double ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, ResidualL2, ResidualL2Action);

    template <typename Archive>
    void save(Archive& ar, unsigned) const;

    template <typename Archive>
    void load(Archive& ar, unsigned);
    HPX_SERIALIZATION_SPLIT_MEMBER();

    void on_migrated();  // Do not rename this is overload member of the base class
};

template <typename ProblemType>
HPXSimulationUnit<ProblemType>::HPXSimulationUnit(const std::string& input_string,
                                                  const uint locality_id,
                                                  const uint submesh_id) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string, locality_id, submesh_id);

    ProblemType::initialize_problem_parameters(input.problem_input);

    input.read_mesh();                         // read mesh meta data
    input.read_bcis();                         // read bc data
    input.read_dbmd(locality_id, submesh_id);  // read distributed boundary meta data

    ProblemType::preprocess_mesh_data(input);

    this->discretization.mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);
    this->communicator        = HPXCommunicator(input.mesh_input.dbmd_data);
    this->stepper             = RKStepper(input.stepper_input);
    this->writer              = Writer<ProblemType>(input.writer_input, locality_id, submesh_id);
    this->parser              = typename ProblemType::ProblemParserType(input, locality_id, submesh_id);

    this->problem_input = input.problem_input;
    this->submesh_model = LoadBalancer::AbstractFactory::create_submesh_model<ProblemType>(
        locality_id, submesh_id, input.load_balancer_input);

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    this->discretization.initialize(input, this->communicator, this->writer);
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::FinishPreprocessor() {
    hpx::future<void> receive_future = this->communicator.ReceivePreprocAll(this->stepper.GetTimestamp());

    this->communicator.SendPreprocAll(this->stepper.GetTimestamp());

    return receive_future.then(
        [this](auto&&) { ProblemType::initialize_data_parallel_post_receive_kernel(this->discretization.mesh); });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Launch() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->discretization.mesh);
    }

    uint n_stages = this->stepper.GetNumStages();

    this->discretization.mesh.CallForEachElement([n_stages](auto& elt) { elt.data.resize(n_stages + 1); });
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Stage() {
    return ProblemType::hpx_stage_kernel(this);
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::SwapStates() {
    this->discretization.mesh.CallForEachElement(
        [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); });

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, this->discretization.mesh);
    }
}

template <typename ProblemType>
double HPXSimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    this->discretization.mesh.CallForEachElement([this, &residual_L2](auto& elt) {
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

    ar& stepper& writer& discretization.mesh& problem_input& parser& communicator& submesh_model;
}

template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::load(Archive& ar, unsigned) {
    ar& stepper& writer& discretization.mesh& problem_input& parser& communicator& submesh_model;

    this->writer.StartLog();

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Arriving on locality " << hpx::get_locality_id() << std::endl;
    }
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::on_migrated() {
    this->discretization.mesh.SetMasters();

    this->discretization.mesh.CallForEachElement([this](auto& elt) {
        using MasterType = typename std::remove_reference<decltype(elt)>::type::ElementMasterType;

        using MasterElementTypes = typename decltype(this->discretization.mesh)::MasterElementTypes;

        MasterType& master_elt =
            std::get<Utilities::index<MasterType, MasterElementTypes>::value>(this->discretization.mesh.GetMasters());

        elt.SetMaster(master_elt);

        elt.Initialize();
    });

    initialize_mesh_interfaces_boundaries<ProblemType, HPXCommunicator>(
        discretization.mesh, problem_input, communicator, writer);
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

    hpx::future<void> Stage() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::StageAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> SwapStates() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::SwapStatesAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<double> ResidualL2() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }
};

#endif