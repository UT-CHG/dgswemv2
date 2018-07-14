#ifndef IHDG_SIMULATION_HPP
#define IHDG_SIMULATION_HPP

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "simulation/writer.hpp"

namespace IHDG {
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

    SparseMatrix<double> delta_local_inv;
    SparseMatrix<double> delta_hat_local;
    DVector<double> rhs_local;

    SparseMatrix<double> delta_global;
    SparseMatrix<double> delta_hat_global;
    DVector<double> rhs_global;

    DMatrix<double> global;

  public:
    Simulation() = default;
    Simulation(const std::string& input_string);

    void Run();
    void ComputeL2Residual();

  private:
    friend ProblemType;
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

    this->mesh.CallForEachElement([this](auto& elt) { elt.data.resize(this->n_stages + 1); });

    ProblemType::initialize_global_problem(this);

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    for (uint step = 1; step <= this->n_steps; ++step) {
        for (uint stage = 0; stage < this->n_stages; ++stage) {
            if (this->parser.ParsingInput()) {
                this->parser.ParseInput(this->stepper, this->mesh);
            }

            ProblemType::initialize_iteration(this);

            uint iter = 0;
            while (true) {
                iter++;

                /* Local Step */
                this->mesh.CallForEachElement(
                    [this](auto& elt) { ProblemType::local_volume_kernel(this->stepper, elt); });

                this->mesh.CallForEachElement(
                    [this](auto& elt) { ProblemType::local_source_kernel(this->stepper, elt); });

                this->mesh.CallForEachInterface(
                    [this](auto& intface) { ProblemType::local_interface_kernel(this->stepper, intface); });

                this->mesh.CallForEachBoundary(
                    [this](auto& bound) { ProblemType::local_boundary_kernel(this->stepper, bound); });

                this->mesh_skeleton.CallForEachEdgeInterface(
                    [this](auto& edge_int) { ProblemType::local_edge_interface_kernel(this->stepper, edge_int); });

                this->mesh_skeleton.CallForEachEdgeBoundary(
                    [this](auto& edge_bound) { ProblemType::local_edge_boundary_kernel(this->stepper, edge_bound); });
                /* Local Step */

                /* Global Step */
                this->mesh_skeleton.CallForEachEdgeInterface(
                    [this](auto& edge_int) { ProblemType::global_edge_interface_kernel(this->stepper, edge_int); });

                this->mesh_skeleton.CallForEachEdgeBoundary(
                    [this](auto& edge_bound) { ProblemType::global_edge_boundary_kernel(this->stepper, edge_bound); });
                /* Global Step */

                bool converged = ProblemType::solve_global_problem(this);

                if (converged) {
                    break;
                }

                if (iter == 100) {
                    break;
                }
            }

            this->mesh.CallForEachElement([this](auto& elt) {
                bool nan_found = ProblemType::scrutinize_solution_kernel(this->stepper, elt);

                if (nan_found)
                    abort();
            });

            ++(this->stepper);
        }

        this->mesh.CallForEachElement([this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); });

        if (this->writer.WritingOutput()) {
            this->writer.WriteOutput(this->stepper, this->mesh);
        }
    }
}

template <typename ProblemType>
void Simulation<ProblemType>::ComputeL2Residual() {
    double residual_L2 = 0;

    this->mesh.CallForEachElement([this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    });

    std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(residual_L2) << std::endl;
}
}

#endif
