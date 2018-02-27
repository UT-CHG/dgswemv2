#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "writer.hpp"

template <typename ProblemType>
class Simulation {
  private:
    Stepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;
    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemInputType problem_input;

  public:
    Simulation(std::string input_string) {
        InputParameters<typename ProblemType::ProblemInputType> input(input_string);

        this->stepper = Stepper(input.rk.nstages, input.rk.order, input.dt, input.T_end);
        this->problem_input = input.problem_input;
        this->writer = Writer<ProblemType>(input);
        this->parser = typename ProblemType::ProblemParserType(input);
        this->mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);

        ProblemType::initialize_problem_parameters(this->problem_input);

        input.ReadMesh();
        mesh.SetMeshName(input.mesh_data.mesh_name);

        std::tuple<> empty_comm;
        initialize_mesh<ProblemType>(this->mesh, input.mesh_data, empty_comm, this->problem_input, this->writer);

        ProblemType::initialize_data_kernel(this->mesh, input.mesh_data, this->problem_input);

        if (this->writer.WritingLog()) {
            this->writer.StartLog();

            this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                      << input.mesh_data.mesh_name << " mesh" << std::endl << std::endl;
        }
    }

    void Run();
    void ComputeL2Residual();
};

template <typename ProblemType>
void Simulation<ProblemType>::Run() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

    auto scrutinize_solution_kernel = [this](auto& elt) {
        ProblemType::scrutinize_solution_kernel(this->stepper, elt);
    };

    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    uint nsteps = (uint)std::ceil(this->stepper.T_end / this->stepper.get_dt());
    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages + 1); };

    this->mesh.CallForEachElement(resize_data_container);

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    for (uint step = 1; step <= nsteps; ++step) {
        if (this->parser.ParsingInput()) {
            this->parser.ParseInput(this->stepper, this->mesh);
        }

        for (uint stage = 0; stage < this->stepper.get_num_stages(); ++stage) {
            this->mesh.CallForEachElement(volume_kernel);

            this->mesh.CallForEachElement(source_kernel);

            this->mesh.CallForEachInterface(interface_kernel);

            this->mesh.CallForEachBoundary(boundary_kernel);

            this->mesh.CallForEachElement(update_kernel);

            this->mesh.CallForEachElement(scrutinize_solution_kernel);

            ProblemType::postprocessor_serial_kernel(this->stepper, this->mesh);

            ++(this->stepper);
        }

        this->mesh.CallForEachElement(swap_states_kernel);

        if (this->writer.WritingOutput()) {
            this->writer.WriteOutput(this->stepper, this->mesh);
        }
    }
}

template <typename ProblemType>
void Simulation<ProblemType>::ComputeL2Residual() {
    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachElement(compute_residual_L2_kernel);

    std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(residual_L2) << std::endl;
}

#endif
