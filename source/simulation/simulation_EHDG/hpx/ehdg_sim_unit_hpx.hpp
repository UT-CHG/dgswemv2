#ifndef EHDG_SIM_UNIT_HPX_HPP
#define EHDG_SIM_UNIT_HPX_HPP

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

  private:
    friend ProblemType;
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

    ProblemType::initialize_global_problem(this);

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
        [this](auto& dbound) { ProblemType::global_distributed_boundary_kernel(this->stepper, dbound); });

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

    /* Global Pre Receive Step */
    this->mesh.CallForEachInterface(
        [this](auto& intface) { ProblemType::global_interface_kernel(this->stepper, intface); });

    this->mesh.CallForEachBoundary([this](auto& bound) { ProblemType::global_boundary_kernel(this->stepper, bound); });

    this->mesh_skeleton.CallForEachEdgeInterface(
        [this](auto& edge_int) { ProblemType::global_edge_interface_kernel(this->stepper, edge_int); });

    this->mesh_skeleton.CallForEachEdgeBoundary(
        [this](auto& edge_bound) { ProblemType::global_edge_boundary_kernel(this->stepper, edge_bound); });
    /* Global Pre Receive Step */

    /* Local Pre Receive Step */
    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::local_volume_kernel(this->stepper, elt); });

    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::local_source_kernel(this->stepper, elt); });

    this->mesh.CallForEachInterface(
        [this](auto& intface) { ProblemType::local_interface_kernel(this->stepper, intface); });

    this->mesh.CallForEachBoundary([this](auto& bound) { ProblemType::local_boundary_kernel(this->stepper, bound); });
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

        /* Global Post Receive Step */
        this->mesh_skeleton.CallForEachEdgeDistributed(
            [this](auto& edge_dbound) { ProblemType::global_edge_distributed_kernel(this->stepper, edge_dbound); });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        this->mesh.CallForEachDistributedBoundary(
            [this](auto& dbound) { ProblemType::local_distributed_boundary_kernel(this->stepper, dbound); });

        this->mesh.CallForEachElement([this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); });

        this->mesh.CallForEachElement([this](auto& elt) {
            bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

            if (nan_found)
                hpx::terminate();
        });
        /* Local Post Receive Step */

        ++(this->stepper);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Step() {
    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); });

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, this->mesh);
    }
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
}

#endif
