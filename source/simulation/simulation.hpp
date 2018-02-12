#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "writer.hpp"

template <typename ProblemType>
class Simulation {
  private:
    InputParameters<typename ProblemType::ProblemInputType> input;

    Stepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemMeshType mesh;

  public:
    Simulation()
        : input(),
          stepper(this->input.rk.nstages, this->input.rk.order, this->input.dt),
          mesh(this->input.polynomial_order) {}
    Simulation(std::string input_string)
        : input(input_string),
          stepper(this->input.rk.nstages, this->input.rk.order, this->input.dt),
          writer(input),
          mesh(this->input.polynomial_order) {
        ProblemType::initialize_problem_parameters(this->input.problem_input);

        this->input.ReadMesh();

        mesh.SetMeshName(this->input.mesh_data.mesh_name);

        std::tuple<> empty_comm;

        if (this->writer.WritingLog()) {
            this->writer.StartLog();

            this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                      << input.mesh_data.mesh_name << " mesh" << std::endl << std::endl;
        }

        initialize_mesh<ProblemType>(
            this->mesh, this->input.mesh_data, empty_comm, this->input.problem_input, this->writer);

        ProblemType::initialize_data_kernel(this->mesh, this->input.mesh_data, this->input.problem_input);
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

    uint nsteps = (uint)std::ceil(this->input.T_end / this->stepper.get_dt());
    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages + 1); };

    this->mesh.CallForEachElement(resize_data_container);

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    ProblemType::parse_source_data(this->stepper, this->mesh, this->input.problem_input);

    for (uint step = 1; step <= nsteps; ++step) {
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
