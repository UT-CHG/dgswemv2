#ifndef SIM_UNIT_OMPI_HPP
#define SIM_UNIT_OMPI_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "communication/ompi_communicator.hpp"
#include "utilities/file_exists.hpp"

#include "simulation/writer.hpp"

template <typename ProblemType>
struct OMPISimulationUnit {
    typename ProblemType::ProblemDiscretizationType discretization;

    OMPICommunicator communicator;
    RKStepper stepper;
    typename ProblemType::ProblemWriterType writer;
    typename ProblemType::ProblemParserType parser;

    typename ProblemType::ProblemInputType problem_input;

    OMPISimulationUnit() = default;
    OMPISimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);
};

template <typename ProblemType>
OMPISimulationUnit<ProblemType>::OMPISimulationUnit(const std::string& input_string,
                                                    const uint locality_id,
                                                    const uint submesh_id) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string, locality_id, submesh_id);

    ProblemType::initialize_problem_parameters(input.problem_input);

    input.read_mesh();                         // read mesh meta data
    input.read_bcis();                         // read bc data
    input.read_dbmd(locality_id, submesh_id);  // read distributed boundary meta data

    ProblemType::preprocess_mesh_data(input);

    this->discretization.mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);
    this->communicator        = OMPICommunicator(input.mesh_input.dbmd_data);
    this->stepper             = RKStepper(input.stepper_input);
    this->writer              = typename ProblemType::ProblemWriterType(input.writer_input, locality_id, submesh_id);
    this->parser              = typename ProblemType::ProblemParserType(input, locality_id, submesh_id);

    this->problem_input = input.problem_input;

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    this->discretization.initialize(input, this->communicator, this->writer);

    this->communicator.InitializeCommunication();
}

#endif
