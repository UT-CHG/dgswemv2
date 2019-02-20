#include "general_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_EHDG/ehdg_swe_problem.hpp"

#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "simulation/serial/simulation.hpp"
#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"

#include "utilities/almost_equal.hpp"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    bool any_error = false;

    typename SWE::EHDG::Problem::ProblemDiscretizationType discretization;
    typename SWE::EHDG::Problem::ProblemGlobalDataType global_data;
    typename SWE::EHDG::Problem::ProblemStepperType stepper;
    typename SWE::EHDG::Problem::ProblemWriterType writer;
    typename SWE::EHDG::Problem::ProblemInputType problem_input;

    std::string input_string = "../../test/files_for_testing/jacobians/input.15";

    InputParameters<typename SWE::EHDG::Problem::ProblemInputType> input(input_string);

    SWE::EHDG::Problem::initialize_problem_parameters(input.problem_input);

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    SWE::EHDG::Problem::preprocess_mesh_data(input);

    discretization.mesh = typename SWE::EHDG::Problem::ProblemMeshType(input.polynomial_order);
    stepper             = typename SWE::EHDG::Problem::ProblemStepperType(input.stepper_input);
    writer              = typename SWE::EHDG::Problem::ProblemWriterType(input.writer_input);
    problem_input       = input.problem_input;

    discretization.initialize(input, writer);

    SWE::EHDG::Problem::preprocessor_serial(discretization, global_data, stepper, problem_input);

    discretization.mesh.CallForEachElement([&](auto& elt) {
        // randomly assign q
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    discretization.mesh.CallForEachInterface(
        [&](auto& intface) { SWE::EHDG::Problem::global_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&](auto& bound) { SWE::EHDG::Problem::global_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        // randomly assign q_hat and delta_q_hat
        HybMatrix<double, SWE::n_variables> q_hat(SWE::n_variables, edge_int.edge_data.get_ndof());
        HybMatrix<double, SWE::n_variables> delta_q_hat(SWE::n_variables, edge_int.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                q_hat(var, dof)       = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
                delta_q_hat(var, dof) = 1.0e-8 * (-1.0 + 2.0 * ((double)rand() / (RAND_MAX)));
            }
        }

        // assign random q_hat
        edge_state.q_hat = q_hat;

        // get current rhs and Jacobian
        edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

        for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
            for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
                submatrix(edge_internal.delta_hat_global,
                          SWE::n_variables * dof_i,
                          SWE::n_variables * dof_j,
                          SWE::n_variables,
                          SWE::n_variables) =
                    reshape<double, SWE::n_variables>(
                        edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
            }

            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
                -edge_int.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }

        // difference
        DynVector<double> diff = edge_internal.rhs_global;
        DynVector<double> diff_est =
            edge_internal.delta_hat_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);

        // assign random q_hat incremented by a random perturbation
        edge_state.q_hat = q_hat + delta_q_hat;

        // find next rhs and diff
        edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            subvector(edge_internal.rhs_global, SWE::n_variables * dof, SWE::n_variables) =
                -edge_int.IntegrationLambda(dof, edge_internal.rhs_global_kernel_at_gp);
        }

        diff -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(diff_est[3 * dof], diff[3 * dof], 1.0e12)) {
                std::cerr << "error in ze" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof] << ' ' << diff[3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 1], diff[3 * dof + 1], 1.0e12)) {
                std::cerr << "error in qx" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 1] << ' ' << diff[3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 2], diff[3 * dof + 2], 1.0e12)) {
                std::cerr << "error in qy" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 2] << ' ' << diff[3 * dof + 2] << std::endl;

                any_error = true;
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        // randomly assign q_hat and delta_q_hat
        HybMatrix<double, SWE::n_variables> q_hat(SWE::n_variables, edge_bound.edge_data.get_ndof());
        HybMatrix<double, SWE::n_variables> delta_q_hat(SWE::n_variables, edge_bound.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                q_hat(var, dof)       = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
                delta_q_hat(var, dof) = 1.0e-8 * (-1.0 + 2.0 * ((double)rand() / (RAND_MAX)));
            }
        }

        // assign random q_hat
        edge_state.q_hat = q_hat;

        // get current rhs and Jacobian
        edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

        for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
            for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
                submatrix(edge_internal.delta_hat_global,
                          SWE::n_variables * dof_i,
                          SWE::n_variables * dof_j,
                          SWE::n_variables,
                          SWE::n_variables) =
                    reshape<double, SWE::n_variables>(
                        edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
            }

            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
                -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }

        // difference
        DynVector<double> diff = edge_internal.rhs_global;
        DynVector<double> diff_est =
            edge_internal.delta_hat_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);

        // assign random q_hat incremented by a random perturbation
        edge_state.q_hat = q_hat + delta_q_hat;

        // find next rhs and diff
        edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            subvector(edge_internal.rhs_global, SWE::n_variables * dof, SWE::n_variables) =
                -edge_bound.IntegrationLambda(dof, edge_internal.rhs_global_kernel_at_gp);
        }

        diff -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(diff_est[3 * dof], diff[3 * dof], 1.0e12)) {
                std::cerr << "error in ze" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof] << ' ' << diff[3 * dof] << std::endl;

                std::cout << diff << std::endl;
                std::cout << diff_est << std::endl;

                abort();

                any_error = true;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 1], diff[3 * dof + 1], 1.0e12)) {
                std::cerr << "error in qx" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 1] << ' ' << diff[3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 2], diff[3 * dof + 2], 1.0e12)) {
                std::cerr << "error in qy" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 2] << ' ' << diff[3 * dof + 2] << std::endl;

                any_error = true;
            }
        }
    });

    if (any_error) {
        return 1;
    }

    return 0;
}
