#ifndef IHDG_SIM_UNIT_OMPI_HPP
#define IHDG_SIM_UNIT_OMPI_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "communication/ompi_communicator.hpp"
#include "utilities/file_exists.hpp"

#include "simulation/writer.hpp"

namespace IHDG {
template <typename ProblemType>
class OMPISimulationUnit {
  private:
    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemMeshSkeletonType mesh_skeleton;

    OMPICommunicator communicator;
    RKStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

  public:
    OMPISimulationUnit() = default;
    OMPISimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);

    void PostReceivePrerocStage();
    void WaitAllPreprocSends();

    void Launch();

    void ExchangeData();
    void PreReceiveStage();
    void PostReceiveStage();
    void WaitAllSends();

    void Step();

    double ResidualL2();
};

template <typename ProblemType>
OMPISimulationUnit<ProblemType>::OMPISimulationUnit(const std::string& input_string,
                                                    const uint locality_id,
                                                    const uint submesh_id) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string, locality_id, submesh_id);

    input.read_mesh();                         // read mesh meta data
    input.read_bcis();                         // read bc data
    input.read_dbmd(locality_id, submesh_id);  // read distributed boundary meta data

    this->mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);

    this->communicator = OMPICommunicator(input.mesh_input.dbmd_data);
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

    initialize_mesh<ProblemType, OMPICommunicator>(this->mesh, input, this->communicator, this->writer);
    initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, this->writer);

    ProblemType::initialize_data_parallel_pre_send_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);

    this->communicator.InitializeCommunication();

    this->communicator.ReceivePreprocAll(this->stepper.GetTimestamp());

    this->communicator.SendPreprocAll(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PostReceivePrerocStage() {
    this->communicator.WaitAllPreprocReceives(this->stepper.GetTimestamp());

    ProblemType::initialize_data_parallel_post_receive_kernel(this->mesh);
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::WaitAllPreprocSends() {
    this->communicator.WaitAllPreprocSends(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::Launch() {
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
void OMPISimulationUnit<ProblemType>::ExchangeData() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.GetTimeAtCurrentStage() << ','
                                  << this->stepper.GetStage() << ')' << std::endl;

        this->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    auto global_distributed_boundary_kernel = [this](auto& dbound) {
        ProblemType::global_distributed_boundary_kernel(this->stepper, dbound);
    };

    this->communicator.ReceiveAll(this->stepper.GetTimestamp());

    this->mesh.CallForEachDistributedBoundary(global_distributed_boundary_kernel);

    this->communicator.SendAll(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PreReceiveStage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

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
        this->writer.GetLogFile() << "Finished work before receive" << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PostReceiveStage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting to wait on receive with timestamp: " << this->stepper.GetTimestamp()
                                  << std::endl;
    }

    this->communicator.WaitAllReceives(this->stepper.GetTimestamp());

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
            MPI_Abort(MPI_COMM_WORLD, 0);
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
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::WaitAllSends() {
    this->communicator.WaitAllSends(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::Step() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, mesh);
    }
}

template <typename ProblemType>
double OMPISimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachElement(compute_residual_L2_kernel);

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}
}

#endif
