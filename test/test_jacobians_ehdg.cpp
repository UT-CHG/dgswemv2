#include "general_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_EHDG/ehdg_swe_problem.hpp"
#include "problem/SWE/discretization_EHDG/kernels_preprocessor/ehdg_swe_kernels_preprocessor.hpp"
#include "problem/SWE/discretization_EHDG/kernels_processor/ehdg_swe_kernels_processor.hpp"
#include "problem/SWE/discretization_EHDG/kernels_postprocessor/ehdg_swe_kernels_postprocessor.hpp"

#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "simulation/serial/simulation.hpp"
#include "simulation/stepper/rk_stepper.hpp"

#include "utilities/almost_equal.hpp"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    std::string input_string = "../../test/files_for_testing/jacobians/input.15";

    InputParameters<typename SWE::EHDG::Problem::ProblemInputType> input(input_string);

    typename SWE::EHDG::Problem::ProblemMeshType mesh;
    typename SWE::EHDG::Problem::ProblemMeshSkeletonType mesh_skeleton;

    RKStepper stepper;
    Writer<SWE::EHDG::Problem> writer;

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    mesh = typename SWE::EHDG::Problem::ProblemMeshType(input.polynomial_order);

    stepper = RKStepper(input.stepper_input);
    writer  = Writer<SWE::EHDG::Problem>(input.writer_input);

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

    mesh_skeleton.CallForEachEdgeInterface([](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                              SWE::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof());
    });

    mesh_skeleton.CallForEachEdgeBoundary([](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof());
    });

    mesh.CallForEachElement([&](auto& elt) {
        // randomly assign q
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];

        for (uint dof = 0; dof < elt.data.get_ndof(); dof++) {
            state.q(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    mesh.CallForEachInterface([&](auto& intface) { SWE::EHDG::Problem::global_interface_kernel(stepper, intface); });

    mesh.CallForEachBoundary([&](auto& bound) { SWE::EHDG::Problem::global_boundary_kernel(stepper, bound); });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); dof++) {
            edge_state.q_hat(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }

        // get current rhs and Jacobian
        edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
                edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp);
        }

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
        // get current rhs and Jacobian

        // randomly assign delta_q_hat and increment
        DynVector<double> delta_q_hat(SWE::n_variables * edge_int.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); dof++) {
            delta_q_hat[3 * dof + SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q_hat[3 * dof + SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q_hat[3 * dof + SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }

        delta_q_hat *= 1.0e-8;  // make it small

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) += delta_q_hat[3 * dof];
            edge_state.q_hat(SWE::Variables::qx, dof) += delta_q_hat[3 * dof + 1];
            edge_state.q_hat(SWE::Variables::qy, dof) += delta_q_hat[3 * dof + 2];
        }
        // randomly assign delta_q_hat and increment

        // find next rhs and substract to find diff
        edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
                edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp);
        }

        edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

        for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) -=
                -edge_int.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }
        // find next rhs and substract to find diff

        // difference estimate
        DynVector<double> diff_est(SWE::n_variables * edge_int.edge_data.get_ndof());

        diff_est = edge_internal.delta_hat_global * delta_q_hat;
        // difference estimate

        // compare
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); dof++) {
            if (!Utilities::almost_equal(diff_est[3 * dof], edge_internal.rhs_global[3 * dof], 1.0e12)) {
                std::cerr << "error in ze" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof] << ' ' << edge_internal.rhs_global[3 * dof]
                          << std::endl;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 1], edge_internal.rhs_global[3 * dof + 1], 1.0e12)) {
                std::cerr << "error in qx" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 1] << ' '
                          << edge_internal.rhs_global[3 * dof + 1] << std::endl;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 2], edge_internal.rhs_global[3 * dof + 2], 1.0e12)) {
                std::cerr << "error in qy" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 2] << ' '
                          << edge_internal.rhs_global[3 * dof + 2] << std::endl;
            }
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); dof++) {
            edge_state.q_hat(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }

        // get current rhs and Jacobian
        edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
                edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
        }

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
        // get current rhs and Jacobian

        // randomly assign delta_q_hat and increment
        DynVector<double> delta_q_hat(SWE::n_variables * edge_bound.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); dof++) {
            delta_q_hat[3 * dof + SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q_hat[3 * dof + SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q_hat[3 * dof + SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }

        delta_q_hat *= 1.0e-8;  // make it small

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) += delta_q_hat[3 * dof];
            edge_state.q_hat(SWE::Variables::qx, dof) += delta_q_hat[3 * dof + 1];
            edge_state.q_hat(SWE::Variables::qy, dof) += delta_q_hat[3 * dof + 2];
        }
        // randomly assign delta_q_hat and increment

        // find next rhs and substract to find diff
        edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
                edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
        }

        edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

        for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) -=
                -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }
        // find next rhs and substract to find diff

        // difference estimate
        DynVector<double> diff_est(SWE::n_variables * edge_bound.edge_data.get_ndof());

        diff_est = edge_internal.delta_hat_global * delta_q_hat;
        // difference estimate

        // compare
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); dof++) {
            if (!Utilities::almost_equal(diff_est[3 * dof], edge_internal.rhs_global[3 * dof], 1.0e10)) {
                std::cerr << "error in ze" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof] << ' ' << edge_internal.rhs_global[3 * dof]
                          << std::endl;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 1], edge_internal.rhs_global[3 * dof + 1], 1.0e10)) {
                std::cerr << "error in qx" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 1] << ' '
                          << edge_internal.rhs_global[3 * dof + 1] << std::endl;
            }

            if (!Utilities::almost_equal(diff_est[3 * dof + 2], edge_internal.rhs_global[3 * dof + 2], 1.0e10)) {
                std::cerr << "error in qy" << std::endl;
                std::cout << std::setprecision(15) << diff_est[3 * dof + 2] << ' '
                          << edge_internal.rhs_global[3 * dof + 2] << std::endl;
            }
        }
    });
}
