#ifndef RKDG_SIM_UNIT_OMPI_HPP
#define RKDG_SIM_UNIT_OMPI_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "communication/ompi_communicator.hpp"
#include "utilities/file_exists.hpp"

#include "simulation/writer.hpp"

namespace RKDG {
template <typename ProblemType>
class OMPISimulationUnit {
  private:
    typename ProblemType::ProblemMeshType mesh;

    OMPICommunicator communicator;
    RKStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

  public:
    OMPISimulationUnit() = default;
    OMPISimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);

    void PostReceivePreprocStage();
    void WaitAllPreprocSends();

    void Launch();

    void ExchangeData();
    void PreReceiveStage();
    void PostReceiveStage();
    void WaitAllSends();

    void ExchangePostprocData();
    void PreReceivePostprocStage();
    void PostReceivePostprocStage();
    void WaitAllPostprocSends();

    void SwapStates();

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

    ProblemType::initialize_data_parallel_pre_send_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);

    this->communicator.InitializeCommunication();

    this->communicator.ReceivePreprocAll(this->stepper.GetTimestamp());

    this->communicator.SendPreprocAll(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PostReceivePreprocStage() {
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

    this->mesh.CallForEachElement([n_stages](auto& elt) { elt.data.resize(n_stages + 1); });
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::ExchangeData() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.GetTimeAtCurrentStage() << ','
                                  << this->stepper.GetStage() << ')' << std::endl;

        this->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    this->communicator.ReceiveAll(this->stepper.GetTimestamp());

    this->mesh.CallForEachDistributedBoundary(
        [this](auto& dbound) { ProblemType::distributed_boundary_send_kernel(this->stepper, dbound); });

    this->communicator.SendAll(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PreReceiveStage() {
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

    this->mesh.CallForEachDistributedBoundary(
        [this](auto& dbound) { ProblemType::distributed_boundary_kernel(this->stepper, dbound); });

    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); });

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::WaitAllSends() {
    this->communicator.WaitAllSends(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::ExchangePostprocData() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging postprocessor data" << std::endl;
    }

    this->communicator.ReceivePostprocAll(this->stepper.GetTimestamp());

    ProblemType::postprocessor_parallel_pre_send_kernel(this->stepper, this->mesh);

    this->communicator.SendPostprocAll(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PreReceivePostprocStage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting postprocessor work before receive" << std::endl;
    }

    ProblemType::postprocessor_parallel_pre_receive_kernel(this->stepper, this->mesh);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished postprocessor work before receive" << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PostReceivePostprocStage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting to wait on postprocessor receive with timestamp: "
                                  << this->stepper.GetTimestamp() << std::endl;
    }

    this->communicator.WaitAllPostprocReceives(this->stepper.GetTimestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting postprocessor work after receive" << std::endl;
    }

    ProblemType::postprocessor_parallel_post_receive_kernel(this->stepper, this->mesh);

    this->mesh.CallForEachElement([this](auto& elt) {
        bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

        if (nan_found)
            MPI_Abort(MPI_COMM_WORLD, 0);
    });

    ++(this->stepper);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished postprocessor work after receive" << std::endl << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::WaitAllPostprocSends() {
    this->communicator.WaitAllPostprocSends(this->stepper.GetTimestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::SwapStates() {
    this->mesh.CallForEachElement([this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); });

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, mesh);
    }
}

template <typename ProblemType>
double OMPISimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    this->mesh.CallForEachElement([this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    });

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}
}

#endif
