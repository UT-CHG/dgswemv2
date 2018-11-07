#ifndef SIM_UNIT_HPX_HPP
#define SIM_UNIT_HPX_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "communication/hpx_communicator.hpp"

#include "simulation/hpx/sim_unit_hpx_base.hpp"
//#include "simulation/hpx/load_balancer/base_model.hpp"
//#include "simulation/hpx/load_balancer/abstract_load_balancer_factory.hpp"

#include "problem/definitions.hpp"
#include "problem/hpx_functions.hpp"

template <typename ProblemType>
struct HPXSimulationUnit : public HPXSimulationUnitBase, hpx::components::managed_component_base<HPXSimulationUnit<ProblemType>> {
    // *** //

    //HPX requires these typedefs to properly disambiguate look ups
    using wrapping_type = typename hpx::components::managed_component_base<HPXSimulationUnit<ProblemType>>::wrapping_type;
    using type_holder = HPXSimulationUnit<ProblemType>;
    using base_type_holder = HPXSimulationUnitBase;

    typename ProblemType::ProblemDiscretizationType discretization;

    HPXCommunicator communicator;
    typename ProblemType::ProblemStepperType stepper;
    typename ProblemType::ProblemWriterType writer;
    typename ProblemType::ProblemParserType parser;

    typename ProblemType::ProblemInputType problem_input;

//    std::unique_ptr<LoadBalancer::SubmeshModel> submesh_model = nullptr;

    HPXSimulationUnit() = default;
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);
    HPXSimulationUnit(HPXSimulationUnit&& rhs) = default;

    HPXSimulationUnit& operator=(HPXSimulationUnit&& rhs) = default;

    hpx::future<void> Preprocessor();

    void Launch();

    hpx::future<void> Step();

    double ResidualL2();

/*    template <typename Archive>
    void save(Archive& ar, unsigned) const;

    template <typename Archive>
    void load(Archive& ar, unsigned);
    HPX_SERIALIZATION_SPLIT_MEMBER();

    void on_migrated();  // Do not rename this is overload member of the base class
*/
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
    this->stepper             = typename ProblemType::ProblemStepperType(input.stepper_input);
    this->writer              = typename ProblemType::ProblemWriterType(input.writer_input, locality_id, submesh_id);
    this->parser              = typename ProblemType::ProblemParserType(input, locality_id, submesh_id);

    this->problem_input = input.problem_input;

/*    this->submesh_model = nullptr; LoadBalancer::AbstractFactory::create_submesh_model<ProblemType>(
                                     locality_id, submesh_id, input.load_balancer_input);*/
    std::cout << "Building sim unit" << '\n';
    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    this->discretization.initialize(input, this->communicator, this->writer);
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Preprocessor() {
    return ProblemType::preprocessor_hpx(this);
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
hpx::future<void> HPXSimulationUnit<ProblemType>::Step() {
    hpx::future<void> step_future = hpx::make_ready_future();

    for (uint stage = 0; stage < this->stepper.GetNumStages(); stage++) {
        step_future = step_future.then([this](auto&& f) {
            f.get();
            if (this->parser.ParsingInput()) {
                this->parser.ParseInput(this->stepper, this->discretization.mesh);
            }
            return ProblemType::stage_hpx(this);
        });
    }

    return step_future.then([this](auto&& f) {
        f.get();
/*        if (this->submesh_model) {
            this->submesh_model->InStep(0, 0);
            }*/

        if (this->writer.WritingOutput()) {
            this->writer.WriteOutput(this->stepper, this->discretization.mesh);
        }
    });
}

template <typename ProblemType>
double HPXSimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    this->discretization.mesh.CallForEachElement(
        [this, &residual_L2](auto& elt) { residual_L2 += ProblemType::compute_residual_L2(this->stepper, elt); });

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}
/*
template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::save(Archive& ar, unsigned) const {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Departing from locality " << hpx::get_locality_id() << std::endl;
    }

    ar& stepper& writer& discretization& problem_input& parser& communicator& submesh_model;
}

template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::load(Archive& ar, unsigned) {
    ar& stepper& writer& discretization& problem_input& parser& communicator& submesh_model;

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
*/

struct HPXEmptySimUnit : HPXSimulationUnitBase, hpx::components::managed_component_base<HPXEmptySimUnit> {
    //HPX requires these typedefs to properly disambiguate look ups
    using wrapping_type = typename hpx::components::managed_component_base<HPXEmptySimUnit>::wrapping_type;
    using type_holder = HPXEmptySimUnit;
    using base_type_holder = HPXSimulationUnitBase;

    hpx::future<void> Preprocessor() { return hpx::make_ready_future(); }

    void Launch() {}

    hpx::future<void> Step() { return hpx::make_ready_future(); }
    double ResidualL2() { return 0.; }
};

using RKDG_SWE_SimUnit = std::conditional<Utilities::is_defined<SWE::RKDG::Problem>::value,
                                          HPXSimulationUnit<SWE::RKDG::Problem>,
                                          HPXEmptySimUnit>::type;
using RKDG_SWE_Server = hpx::components::managed_component<RKDG_SWE_SimUnit>;
HPX_REGISTER_DERIVED_COMPONENT_FACTORY(RKDG_SWE_Server, RKDG_SWE_SimUnit, "HPXSimulationUnitBase");

using EHDG_SWE_SimUnit = std::conditional<Utilities::is_defined<SWE::EHDG::Problem>::value,
                                          HPXSimulationUnit<SWE::EHDG::Problem>,
                                          HPXEmptySimUnit>::type;
using EHDG_SWE_Server = hpx::components::managed_component<EHDG_SWE_SimUnit>;
HPX_REGISTER_DERIVED_COMPONENT_FACTORY(EHDG_SWE_Server, EHDG_SWE_SimUnit, "HPXSimulationUnitBase");

#endif