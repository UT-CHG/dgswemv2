#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "preprocessor/input_parameters.hpp"
#include "simulation/writer.hpp"

template <typename ProblemType>
class Simulation {
  private:
    uint n_steps;
    uint n_stages;

    typename ProblemType::ProblemDiscretizationType discretization;

    RKStepper stepper;
    typename ProblemType::ProblemWriterType writer;
    typename ProblemType::ProblemParserType parser;

    typename ProblemType::ProblemInputType problem_input;

  public:
    Simulation() = default;
    Simulation(const std::string& input_string);

    void Run();
    void ComputeL2Residual();
};

template <typename ProblemType>
Simulation<ProblemType>::Simulation(const std::string& input_string) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string);

    ProblemType::initialize_problem_parameters(input.problem_input);

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    ProblemType::preprocess_mesh_data(input);

    this->n_steps  = (uint)std::ceil(input.stepper_input.run_time / input.stepper_input.dt);
    this->n_stages = input.stepper_input.nstages;

    this->discretization.mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);
    this->stepper             = RKStepper(input.stepper_input);
    this->writer              = typename ProblemType::ProblemWriterType(input.writer_input);
    this->parser              = typename ProblemType::ProblemParserType(input);

    this->problem_input = input.problem_input;

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    this->discretization.initialize(input, this->writer);
}

template <typename ProblemType>
void Simulation<ProblemType>::Run() {
    ProblemType::preprocessor_serial(this->discretization, this->problem_input);

    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    this->discretization.mesh.CallForEachElement([this](auto& elt) { elt.data.resize(this->n_stages + 1); });

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->discretization.mesh);
    }

    for (uint step = 1; step <= this->n_steps; ++step) {
        for (uint stage = 0; stage < this->n_stages; ++stage) {
            if (this->parser.ParsingInput()) {
                this->parser.ParseInput(this->stepper, this->discretization.mesh);
            }

            ProblemType::stage_serial(this->stepper, this->discretization);

            ++(this->stepper);
        }

        this->discretization.mesh.CallForEachElement(
            [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); });

        if (this->writer.WritingOutput()) {
            this->writer.WriteOutput(this->stepper, this->discretization.mesh);
        }
    }
}

template <typename ProblemType>
void Simulation<ProblemType>::ComputeL2Residual() {
    double residual_L2 = 0;

    this->discretization.mesh.CallForEachElement(
        [this, &residual_L2](auto& elt) { residual_L2 += ProblemType::compute_residual_L2(this->stepper, elt); });

    std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(residual_L2) << std::endl;
}

#endif
