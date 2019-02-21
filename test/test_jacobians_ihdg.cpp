#include "general_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_IHDG/ihdg_swe_problem.hpp"

#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "simulation/serial/simulation.hpp"
#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"

#include "utilities/almost_equal.hpp"

std::vector<std::string> var_name = {"ze", "qx", "qy"};

int main(int argc, char* argv[]) {
    srand(time(NULL));

    bool any_error = false;

    typename SWE::IHDG::Problem::ProblemDiscretizationType discretization;
    typename SWE::IHDG::Problem::ProblemGlobalDataType global_data;
    typename SWE::IHDG::Problem::ProblemStepperType stepper;
    typename SWE::IHDG::Problem::ProblemWriterType writer;
    typename SWE::IHDG::Problem::ProblemInputType problem_input;

    std::string input_string = "../../test/files_for_testing/jacobians/input.15";

    InputParameters<typename SWE::IHDG::Problem::ProblemInputType> input(input_string);

    SWE::IHDG::Problem::initialize_problem_parameters(input.problem_input);

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    SWE::IHDG::Problem::preprocess_mesh_data(input);

    discretization.mesh = typename SWE::IHDG::Problem::ProblemMeshType(input.polynomial_order);
    stepper             = typename SWE::IHDG::Problem::ProblemStepperType(input.stepper_input);
    writer              = typename SWE::IHDG::Problem::ProblemWriterType(input.writer_input);
    problem_input       = input.problem_input;

    discretization.initialize(input, writer);

    SWE::IHDG::Problem::preprocessor_serial(discretization, global_data, stepper, problem_input);

    discretization.mesh.CallForEachElement([&](auto& elt) {
        elt.data.resize(2);

        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];

        // randomly assign q
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                state.q(var, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                edge_state.q_hat(var, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                edge_state.q_hat(var, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            }
        }
    });

    SWE::IHDG::Problem::init_iteration(stepper, discretization);

    // do one pass to compute all jacobians and rhs
    discretization.mesh.CallForEachElement([&](auto& elt) {
        SWE::IHDG::Problem::local_volume_kernel(stepper, elt);
        SWE::IHDG::Problem::local_source_kernel(stepper, elt);
    });

    discretization.mesh.CallForEachInterface(
        [&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    // do one pass to compute all jacobians and rhs

    // containers to store data
    std::vector<DynVector<double>> delta_diff(discretization.mesh.GetNumberElements());
    std::vector<DynVector<double>> delta_diff_est(discretization.mesh.GetNumberElements());

    discretization.mesh.CallForEachElement([&](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage + 1];

        // randomly assign delta_q and increment
        HybMatrix<double, SWE::n_variables> delta_q(SWE::n_variables, elt.data.get_ndof());

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                delta_q(var, dof) = 1.0e-8 * (-1.0 + 2.0 * ((double)rand() / (RAND_MAX)));
            }
        }

        // assign random q incremented by a random perturbation
        state.q += delta_q;

        // store rhs_local for future comparison
        delta_diff[elt.GetID()] = elt.data.internal.rhs_local;

        // difference estimate
        delta_diff_est[elt.GetID()] =
            elt.data.internal.delta_local * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q);

        // store for global jacobian checks
        elt.data.state[0].q = delta_q;
    });

    // containers to store data
    std::vector<DynVector<double>> delta_hat_diff(discretization.mesh_skeleton.GetNumberEdgeBoundaries() +
                                                  discretization.mesh_skeleton.GetNumberEdgeInterfaces());
    std::vector<DynVector<double>> delta_hat_diff_est(discretization.mesh_skeleton.GetNumberEdgeBoundaries() +
                                                      discretization.mesh_skeleton.GetNumberEdgeInterfaces());
    // containers to store data

    discretization.mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& state_in = edge_int.interface.data_in.state[0];
        auto& state_ex = edge_int.interface.data_ex.state[0];

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint edge_ID = edge_internal.global_dof_indx[0] / edge_int.edge_data.get_ndof() / SWE::n_variables;

        // store rhs_global for future comparison
        delta_hat_diff[edge_ID] = edge_internal.rhs_global;

        // difference estimate
        delta_hat_diff_est[edge_ID] =
            boundary_in.delta_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(state_in.q) +
            boundary_ex.delta_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(state_ex.q);
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& state = edge_bound.boundary.data.state[0];

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint edge_ID = edge_internal.global_dof_indx[0] / edge_bound.edge_data.get_ndof() / SWE::n_variables;

        // store rhs_global for future comparison
        delta_hat_diff[edge_ID] = edge_internal.rhs_global;

        // difference estimate
        delta_hat_diff_est[edge_ID] =
            boundary.delta_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(state.q);
    });

    // do one pass to compute next rhs
    discretization.mesh.CallForEachElement([&](auto& elt) {
        SWE::IHDG::Problem::local_volume_kernel(stepper, elt);
        SWE::IHDG::Problem::local_source_kernel(stepper, elt);
    });

    discretization.mesh.CallForEachInterface(
        [&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh.CallForEachElement([&](auto& elt) {
        // subtract to find diff true
        delta_diff[elt.GetID()] -= elt.data.internal.rhs_local;

        // compare
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                if (!Utilities::almost_equal(delta_diff_est[elt.GetID()][SWE::n_variables * dof + var],
                                             delta_diff[elt.GetID()][SWE::n_variables * dof + var],
                                             1.0e12)) {
                    std::cerr << "error del_local in " << var_name[var] << std::endl;
                    std::cout << std::setprecision(15) << delta_diff_est[elt.GetID()][SWE::n_variables * dof + var]
                              << ' ' << delta_diff[elt.GetID()][SWE::n_variables * dof + var] << std::endl;

                    any_error = true;
                }
            }
        }

        // store rhs_local for future comparison
        delta_diff[elt.GetID()] = elt.data.internal.rhs_local;

        // difference estimate zero out for future
        set_constant(elt.data.internal.rhs_local, 0.0);
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint edge_ID = edge_internal.global_dof_indx[0] / edge_int.edge_data.get_ndof() / SWE::n_variables;

        // subtract to find diff true
        delta_hat_diff[edge_ID] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                if (!Utilities::almost_equal(delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var],
                                             delta_hat_diff[edge_ID][SWE::n_variables * dof + var],
                                             1.0e12)) {
                    std::cerr << "error edge int del_global in " << var_name[var] << std::endl;
                    std::cout << std::setprecision(15) << delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var]
                              << ' ' << delta_hat_diff[edge_ID][SWE::n_variables * dof + var] << std::endl;

                    any_error = true;
                }
            }
        }

        // randomly assign delta_q_hat and increment
        HybMatrix<double, SWE::n_variables> delta_q_hat(SWE::n_variables, edge_int.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                delta_q_hat(var, dof) = 1.0e-8 * (-1.0 + 2.0 * ((double)rand() / (RAND_MAX)));
            }
        }

        edge_state.q_hat += delta_q_hat;

        // difference estimate
        delta_hat_diff[edge_ID] = edge_internal.rhs_global;

        delta_hat_diff_est[edge_ID] =
            edge_internal.delta_hat_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);

        internal_in.rhs_local +=
            boundary_in.delta_hat_local * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);
        internal_ex.rhs_local +=
            boundary_ex.delta_hat_local * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint edge_ID = edge_internal.global_dof_indx[0] / edge_bound.edge_data.get_ndof() / SWE::n_variables;

        // subtract to find diff true
        delta_hat_diff[edge_ID] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                if (!Utilities::almost_equal(delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var],
                                             delta_hat_diff[edge_ID][SWE::n_variables * dof + var],
                                             1.0e12)) {
                    std::cerr << "error edge bound del_global in " << var_name[var] << std::endl;
                    std::cout << std::setprecision(15) << delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var]
                              << ' ' << delta_hat_diff[edge_ID][SWE::n_variables * dof + var] << std::endl;

                    any_error = true;
                }
            }
        }

        // randomly assign delta_q_hat and increment
        HybMatrix<double, SWE::n_variables> delta_q_hat(SWE::n_variables, edge_bound.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                delta_q_hat(var, dof) = 1.0e-8 * (-1.0 + 2.0 * ((double)rand() / (RAND_MAX)));
            }
        }

        edge_state.q_hat += delta_q_hat;

        delta_hat_diff[edge_ID] = edge_internal.rhs_global;

        delta_hat_diff_est[edge_ID] =
            edge_internal.delta_hat_global * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);

        internal.rhs_local +=
            boundary.delta_hat_local * flatten<double, SWE::n_variables, SO::ColumnMajor>(delta_q_hat);
    });

    discretization.mesh.CallForEachElement([&](auto& elt) {
        // difference estimate zero out for future
        delta_diff_est[elt.GetID()] = elt.data.internal.rhs_local;
    });

    // do one pass to compute next rhs
    discretization.mesh.CallForEachElement([&](auto& elt) { SWE::IHDG::Problem::local_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    // do one pass to compute next rhs

    discretization.mesh.CallForEachElement([&](auto& elt) {
        // subtract to find diff true
        delta_diff[elt.GetID()] -= elt.data.internal.rhs_local;

        // compare
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                if (!Utilities::almost_equal(delta_diff_est[elt.GetID()][SWE::n_variables * dof + var],
                                             delta_diff[elt.GetID()][SWE::n_variables * dof + var],
                                             1.0e12)) {
                    std::cerr << "error del_hat_local in ze" << std::endl;
                    std::cout << std::setprecision(15) << delta_diff_est[elt.GetID()][SWE::n_variables * dof + var]
                              << ' ' << delta_diff[elt.GetID()][SWE::n_variables * dof + var] << std::endl;

                    any_error = true;
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        uint edge_ID = edge_internal.global_dof_indx[0] / edge_int.edge_data.get_ndof() / SWE::n_variables;

        // subtract to find diff true
        delta_hat_diff[edge_ID] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                if (!Utilities::almost_equal(delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var],
                                             delta_hat_diff[edge_ID][SWE::n_variables * dof + var],
                                             1.0e12)) {
                    std::cerr << "error edge int del_hat_global in ze" << std::endl;
                    std::cout << std::setprecision(15) << delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var]
                              << ' ' << delta_hat_diff[edge_ID][SWE::n_variables * dof + var] << std::endl;

                    any_error = true;
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        uint edge_ID = edge_internal.global_dof_indx[0] / edge_bound.edge_data.get_ndof() / SWE::n_variables;

        // subtract to find diff true
        delta_hat_diff[edge_ID] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                if (!Utilities::almost_equal(delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var],
                                             delta_hat_diff[edge_ID][SWE::n_variables * dof + var],
                                             1.0e12)) {
                    std::cerr << "error edge bound del_hat_global in ze" << std::endl;
                    std::cout << std::setprecision(15) << delta_hat_diff_est[edge_ID][SWE::n_variables * dof + var]
                              << ' ' << delta_hat_diff[edge_ID][SWE::n_variables * dof + var] << std::endl;

                    any_error = true;
                }
            }
        }
    });

    if (any_error) {
        return 1;
    }

    return 0;
}