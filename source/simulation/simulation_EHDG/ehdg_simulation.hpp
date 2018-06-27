#ifndef EHDG_SIMULATION_HPP
#define EHDG_SIMULATION_HPP

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "simulation/writer.hpp"

namespace EHDG {
template <typename ProblemType>
class Simulation {
  private:
    uint n_steps;
    uint n_stages;

    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemMeshSkeletonType mesh_skeleton;

    RKStepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;

  public:
    Simulation() = default;
    Simulation(const std::string& input_string);

    void Run();
    void ComputeL2Residual();
};

template <typename ProblemType>
Simulation<ProblemType>::Simulation(const std::string& input_string) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string);

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    this->n_steps  = (uint)std::ceil(input.stepper_input.run_time / input.stepper_input.dt);
    this->n_stages = input.stepper_input.nstages;

    this->mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);

    this->stepper = RKStepper(input.stepper_input);
    this->writer  = Writer<ProblemType>(input.writer_input);
    this->parser  = typename ProblemType::ProblemParserType(input);

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                                  << std::endl;
    }

    ProblemType::initialize_problem_parameters(input.problem_input);

    ProblemType::preprocess_mesh_data(input);

    std::tuple<> empty_comm;

    initialize_mesh<ProblemType>(this->mesh, input, empty_comm, this->writer);
    initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, this->writer);

    ProblemType::initialize_data_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);
}

template <typename ProblemType>
void Simulation<ProblemType>::Run() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    auto global_interface_kernel = [this](auto& intface) {
        ProblemType::global_interface_kernel(this->stepper, intface);
    };

    auto global_boundary_kernel = [this](auto& bound) { ProblemType::global_boundary_kernel(this->stepper, bound); };

    auto global_edge_internal_kernel = [this](auto& edge_int) {
        ProblemType::global_edge_internal_kernel(this->stepper, edge_int);
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

    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

    auto scrutinize_solution_kernel = [this](auto& elt) {
        bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

        if (nan_found)
            abort();
    };

    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    auto resize_data_container = [this](auto& elt) { elt.data.resize(this->n_stages + 1); };

    this->mesh.CallForEachElement(resize_data_container);

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    for (uint step = 1; step <= this->n_steps; ++step) {
        for (uint stage = 0; stage < this->n_stages; ++stage) {
            if (this->parser.ParsingInput()) {
                this->parser.ParseInput(this->stepper, this->mesh);
            }
            /* Global Step */
            this->mesh.CallForEachInterface(global_interface_kernel);

            // this->mesh.CallForEachBoundary(global_boundary_kernel);

            this->mesh_skeleton.CallForEachEdgeInternal(global_edge_internal_kernel);

            // this->mesh_skeleton.CallForEachEdgeBoundary(global_edge_boundary_kernel);
            /* Global Step */

            /* Local Step */
            this->mesh.CallForEachElement(local_volume_kernel);

            this->mesh.CallForEachElement(local_source_kernel);

            this->mesh.CallForEachInterface(local_interface_kernel);

            this->mesh.CallForEachBoundary(local_boundary_kernel);

            this->mesh.CallForEachElement(update_kernel);

            this->mesh.CallForEachElement(scrutinize_solution_kernel);
            /* Local Step */

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
}

#endif
