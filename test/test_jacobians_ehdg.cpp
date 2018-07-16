#include "general_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_EHDG/ehdg_swe_problem.hpp"
#include "problem/SWE/discretization_EHDG/kernels_preprocessor/ehdg_swe_kernels_preprocessor.hpp"
#include "problem/SWE/discretization_EHDG/kernels_processor/ehdg_swe_kernels_processor.hpp"
#include "problem/SWE/discretization_EHDG/kernels_postprocessor/ehdg_swe_kernels_postprocessor.hpp"

#include "simulation/simulation_EHDG/serial/ehdg_simulation.hpp"
#include "simulation/stepper/rk_stepper.hpp"

int main(int argc, char* argv[]) {
    std::string input_string = "../../test/files_for_testing/jacobians/input.15";

    InputParameters<typename SWE::EHDG::Problem::ProblemInputType> input(input_string);

    typename SWE::EHDG::Problem::ProblemMeshType mesh;
    typename SWE::EHDG::Problem::ProblemMeshSkeletonType mesh_skeleton;

    Writer<SWE::EHDG::Problem> writer;

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    mesh = typename SWE::EHDG::Problem::ProblemMeshType(input.polynomial_order);

    writer = Writer<SWE::EHDG::Problem>(input.writer_input);

    if (writer.WritingLog()) {
        writer.StartLog();

        writer.GetLogFile() << "Starting check jacobian with p=" << input.polynomial_order << " for "
                            << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                            << std::endl;
    }

    SWE::EHDG::Problem::initialize_problem_parameters(input.problem_input);

    SWE::EHDG::Problem::preprocess_mesh_data(input);

    std::tuple<> empty_comm;

    initialize_mesh<SWE::EHDG::Problem>(mesh, input, empty_comm, writer);
    initialize_mesh_skeleton<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);

    SWE::EHDG::Problem::initialize_data_kernel(mesh, input.mesh_input.mesh_data, input.problem_input);
}
