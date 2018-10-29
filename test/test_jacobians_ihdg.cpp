#include "general_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_IHDG/ihdg_swe_problem.hpp"
#include "problem/SWE/discretization_IHDG/kernels_preprocessor/ihdg_swe_kernels_preprocessor.hpp"
#include "problem/SWE/discretization_IHDG/kernels_processor/ihdg_swe_kernels_processor.hpp"
#include "problem/SWE/discretization_IHDG/kernels_postprocessor/ihdg_swe_kernels_postprocessor.hpp"

#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "simulation/serial/simulation.hpp"
#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"

#include "utilities/almost_equal.hpp"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    bool any_error = false;

    std::string input_string = "../../test/files_for_testing/jacobians/input.15";

    InputParameters<typename SWE::IHDG::Problem::ProblemInputType> input(input_string);

    typename SWE::IHDG::Problem::ProblemMeshType mesh;
    typename SWE::IHDG::Problem::ProblemMeshSkeletonType mesh_skeleton;

    ESSPRKStepper stepper;
    Writer<SWE::IHDG::Problem> writer;

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    mesh = typename SWE::IHDG::Problem::ProblemMeshType(input.polynomial_order);

    stepper = ESSPRKStepper(input.stepper_input);
    writer  = Writer<SWE::IHDG::Problem>(input.writer_input);

    if (writer.WritingLog()) {
        writer.StartLog();

        writer.GetLogFile() << "Starting check jacobian with p=" << input.polynomial_order << " for "
                            << input.mesh_input.mesh_data.mesh_name << " mesh" << std::endl
                            << std::endl;
    }

    SWE::IHDG::Problem::initialize_problem_parameters(input.problem_input);

    SWE::IHDG::Problem::preprocess_mesh_data(input);

    std::tuple<> empty_comm;

    initialize_mesh<SWE::IHDG::Problem>(mesh, input, empty_comm, writer);
    initialize_mesh_skeleton<SWE::IHDG::Problem>(mesh, mesh_skeleton, writer);

    SWE::IHDG::Problem::initialize_data_serial(mesh, input.problem_input);

    mesh.CallForEachElement([](auto& elt) { elt.data.resize(2); });

    uint local_dof_offset  = 0;
    uint global_dof_offset = 0;

    mesh.CallForEachElement([&](auto& elt) {
        auto& internal = elt.data.internal;

        // Set offsets for global matrix construction
        internal.local_dof_offset = local_dof_offset;
        local_dof_offset += elt.data.get_ndof() * SWE::n_variables;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].local_dof_offset = internal.local_dof_offset;
        }

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;
        global_dof_offset += edge_int.edge_data.get_ndof() * SWE::n_variables;

        boundary_in.global_dof_offset = global_dof_offset;
        boundary_ex.global_dof_offset = global_dof_offset;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                              SWE::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_hat_local containers
        boundary_in.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_in.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());
        boundary_ex.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_ex.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_global containers
        boundary_in.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_in.get_ndof());
        boundary_ex.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_ex.get_ndof());
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;
        global_dof_offset += edge_bound.edge_data.get_ndof() * SWE::n_variables;

        boundary.global_dof_offset = global_dof_offset;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_hat_local container
        boundary.delta_hat_local.resize(SWE::n_variables * edge_bound.boundary.data.get_ndof(),
                                        SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_global container
        boundary.delta_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                     SWE::n_variables * edge_bound.boundary.data.get_ndof());
    });

    mesh.CallForEachElement([&](auto& elt) {
        // randomly assign q
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage + 1];

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qx, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat(SWE::Variables::qy, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    // do one pass to compute all jacobians and rhs
    mesh.CallForEachElement([&](auto& elt) { SWE::IHDG::Problem::local_volume_kernel(stepper, elt); });

    mesh.CallForEachInterface([&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    mesh.CallForEachBoundary([&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    // do one pass to compute all jacobians and rhs

    // containers to store data
    DynVector<DynVector<double>> rhs_local_prev(mesh.GetNumberElements());
    DynVector<DynVector<double>> delta_local_diff_est(mesh.GetNumberElements());
    DynVector<DynVector<double>> delta_q_store(mesh.GetNumberElements());
    // containers to store data

    mesh.CallForEachElement([&](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage + 1];

        // randomly assign delta_q and increment
        DynVector<double> delta_q(SWE::n_variables * elt.data.get_ndof());

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            delta_q[3 * dof + SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q[3 * dof + SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q[3 * dof + SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }

        delta_q *= 1.0e-8;  // make it small

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q(SWE::Variables::ze, dof) += delta_q[3 * dof];
            state.q(SWE::Variables::qx, dof) += delta_q[3 * dof + 1];
            state.q(SWE::Variables::qy, dof) += delta_q[3 * dof + 2];
        }
        // randomly assign delta_q and increment

        // store for global jacobian checks
        delta_q_store[elt.GetID()] = delta_q;

        // difference estimate
        delta_local_diff_est[elt.GetID()] = elt.data.internal.delta_local * delta_q;

        // store rhs_local for future comparison
        rhs_local_prev[elt.GetID()] = elt.data.internal.rhs_local;
    });

    // containers to store data
    DynVector<DynVector<double>> rhs_global_prev(mesh_skeleton.GetNumberEdgeBoundaries() +
                                                 mesh_skeleton.GetNumberEdgeInterfaces());
    DynVector<DynVector<double>> delta_global_diff_est(mesh_skeleton.GetNumberEdgeBoundaries() +
                                                       mesh_skeleton.GetNumberEdgeInterfaces());
    // containers to store data

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // difference estimate
        delta_global_diff_est[edge_int.GetID()] =
            boundary_in.delta_global * delta_q_store[boundary_in.local_dof_offset / 9];
        delta_global_diff_est[edge_int.GetID()] +=
            boundary_ex.delta_global * delta_q_store[boundary_ex.local_dof_offset / 9];

        // store rhs_global for future comparison
        rhs_global_prev[edge_int.GetID()] = edge_internal.rhs_global;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // difference estimate
        delta_global_diff_est[edge_bound.GetID()] =
            boundary.delta_global * delta_q_store[boundary.local_dof_offset / 9];

        // store rhs_global for future comparison
        rhs_global_prev[edge_bound.GetID()] = edge_internal.rhs_global;
    });

    // do one pass to compute next rhs
    mesh.CallForEachElement([&](auto& elt) { SWE::IHDG::Problem::local_volume_kernel(stepper, elt); });

    mesh.CallForEachInterface([&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    mesh.CallForEachBoundary([&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    // do one pass to compute next rhs

    mesh.CallForEachElement([&](auto& elt) {
        // subtract to find diff true
        rhs_local_prev[elt.GetID()] -= elt.data.internal.rhs_local;

        // compare
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(
                    delta_local_diff_est[elt.GetID()][3 * dof], rhs_local_prev[elt.GetID()][3 * dof], 1.0e12)) {
                std::cerr << "error del_local in ze" << std::endl;
                std::cout << std::setprecision(15) << delta_local_diff_est[elt.GetID()][3 * dof] << ' '
                          << rhs_local_prev[elt.GetID()][3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(
                    delta_local_diff_est[elt.GetID()][3 * dof + 1], rhs_local_prev[elt.GetID()][3 * dof + 1], 1.0e12)) {
                std::cerr << "error del_local in qx" << std::endl;
                std::cout << std::setprecision(15) << delta_local_diff_est[elt.GetID()][3 * dof + 1] << ' '
                          << rhs_local_prev[elt.GetID()][3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(
                    delta_local_diff_est[elt.GetID()][3 * dof + 2], rhs_local_prev[elt.GetID()][3 * dof + 2], 1.0e12)) {
                std::cerr << "error del_local in qy" << std::endl;
                std::cout << std::setprecision(15) << delta_local_diff_est[elt.GetID()][3 * dof + 2] << ' '
                          << rhs_local_prev[elt.GetID()][3 * dof + 2] << std::endl;

                any_error = true;
            }
        }

        // difference estimate zero out for future
        set_constant(delta_local_diff_est[elt.GetID()], 0.0);

        // store rhs_local for future comparison
        rhs_local_prev[elt.GetID()] = elt.data.internal.rhs_local;
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // subtract to find diff true
        rhs_global_prev[edge_int.GetID()] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(delta_global_diff_est[edge_int.GetID()][3 * dof],
                                         rhs_global_prev[edge_int.GetID()][3 * dof],
                                         1.0e12)) {
                std::cerr << "error edge int del_global in ze" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_int.GetID()][3 * dof] << ' '
                          << rhs_global_prev[edge_int.GetID()][3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_int.GetID()][3 * dof + 1],
                                         rhs_global_prev[edge_int.GetID()][3 * dof + 1],
                                         1.0e12)) {
                std::cerr << "error edge int del_global in qx" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_int.GetID()][3 * dof + 1] << ' '
                          << rhs_global_prev[edge_int.GetID()][3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_int.GetID()][3 * dof + 2],
                                         rhs_global_prev[edge_int.GetID()][3 * dof + 2],
                                         1.0e12)) {
                std::cerr << "error edge int del_global in qy" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_int.GetID()][3 * dof + 2] << ' '
                          << rhs_global_prev[edge_int.GetID()][3 * dof + 2] << std::endl;

                any_error = true;
            }
        }

        // randomly assign delta_q_hat and increment
        DynVector<double> delta_q_hat(SWE::n_variables * edge_int.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
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

        // difference estimate
        rhs_global_prev[edge_int.GetID()] = edge_internal.rhs_global;

        delta_global_diff_est[edge_int.GetID()] = edge_internal.delta_hat_global * delta_q_hat;

        delta_local_diff_est[boundary_in.local_dof_offset / 9] += boundary_in.delta_hat_local * delta_q_hat;
        delta_local_diff_est[boundary_ex.local_dof_offset / 9] += boundary_ex.delta_hat_local * delta_q_hat;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // subtract to find diff true
        rhs_global_prev[edge_bound.GetID()] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(delta_global_diff_est[edge_bound.GetID()][3 * dof],
                                         rhs_global_prev[edge_bound.GetID()][3 * dof],
                                         1.0e12)) {
                std::cerr << "error edge bound del_global in ze" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_bound.GetID()][3 * dof] << ' '
                          << rhs_global_prev[edge_bound.GetID()][3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_bound.GetID()][3 * dof + 1],
                                         rhs_global_prev[edge_bound.GetID()][3 * dof + 1],
                                         1.0e12)) {
                std::cerr << "error edge bound del_global in qx" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_bound.GetID()][3 * dof + 1] << ' '
                          << rhs_global_prev[edge_bound.GetID()][3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_bound.GetID()][3 * dof + 2],
                                         rhs_global_prev[edge_bound.GetID()][3 * dof + 2],
                                         1.0e12)) {
                std::cerr << "error edge bound del_global in qy" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_bound.GetID()][3 * dof + 2] << ' '
                          << rhs_global_prev[edge_bound.GetID()][3 * dof + 2] << std::endl;

                any_error = true;
            }
        }

        // randomly assign delta_q_hat and increment
        DynVector<double> delta_q_hat(SWE::n_variables * edge_bound.edge_data.get_ndof());

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
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

        rhs_global_prev[edge_bound.GetID()] = edge_internal.rhs_global;

        delta_global_diff_est[edge_bound.GetID()] = edge_internal.delta_hat_global * delta_q_hat;

        // TEMP FIX HERE! boundary.local_dof_offset / 9
        delta_local_diff_est[boundary.local_dof_offset / 9] += boundary.delta_hat_local * delta_q_hat;
    });

    // do one pass to compute next rhs
    mesh.CallForEachElement([&](auto& elt) { SWE::IHDG::Problem::local_volume_kernel(stepper, elt); });

    mesh.CallForEachInterface([&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    mesh.CallForEachBoundary([&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    // do one pass to compute next rhs

    mesh.CallForEachElement([&](auto& elt) {
        // subtract to find diff true
        rhs_local_prev[elt.GetID()] -= elt.data.internal.rhs_local;

        // compare
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(
                    delta_local_diff_est[elt.GetID()][3 * dof], rhs_local_prev[elt.GetID()][3 * dof], 1.0e12)) {
                std::cerr << "error del_hat_local in ze" << std::endl;
                std::cout << std::setprecision(15) << delta_local_diff_est[elt.GetID()][3 * dof] << ' '
                          << rhs_local_prev[elt.GetID()][3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(
                    delta_local_diff_est[elt.GetID()][3 * dof + 1], rhs_local_prev[elt.GetID()][3 * dof + 1], 1.0e12)) {
                std::cerr << "error del_hat_local in qx" << std::endl;
                std::cout << std::setprecision(15) << delta_local_diff_est[elt.GetID()][3 * dof + 1] << ' '
                          << rhs_local_prev[elt.GetID()][3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(
                    delta_local_diff_est[elt.GetID()][3 * dof + 2], rhs_local_prev[elt.GetID()][3 * dof + 2], 1.0e12)) {
                std::cerr << "error del_hat_local in qy" << std::endl;
                std::cout << std::setprecision(15) << delta_local_diff_est[elt.GetID()][3 * dof + 2] << ' '
                          << rhs_local_prev[elt.GetID()][3 * dof + 2] << std::endl;

                any_error = true;
            }
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        // subtract to find diff true
        rhs_global_prev[edge_int.GetID()] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(delta_global_diff_est[edge_int.GetID()][3 * dof],
                                         rhs_global_prev[edge_int.GetID()][3 * dof],
                                         1.0e12)) {
                std::cerr << "error edge int del_hat_global in ze" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_int.GetID()][3 * dof] << ' '
                          << rhs_global_prev[edge_int.GetID()][3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_int.GetID()][3 * dof + 1],
                                         rhs_global_prev[edge_int.GetID()][3 * dof + 1],
                                         1.0e12)) {
                std::cerr << "error edge int del_hat_global in qx" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_int.GetID()][3 * dof + 1] << ' '
                          << rhs_global_prev[edge_int.GetID()][3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_int.GetID()][3 * dof + 2],
                                         rhs_global_prev[edge_int.GetID()][3 * dof + 2],
                                         1.0e12)) {
                std::cerr << "error edge int del_hat_global in qy" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_int.GetID()][3 * dof + 2] << ' '
                          << rhs_global_prev[edge_int.GetID()][3 * dof + 2] << std::endl;

                any_error = true;
            }
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        // subtract to find diff true
        rhs_global_prev[edge_bound.GetID()] -= edge_internal.rhs_global;

        // compare
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            if (!Utilities::almost_equal(delta_global_diff_est[edge_bound.GetID()][3 * dof],
                                         rhs_global_prev[edge_bound.GetID()][3 * dof],
                                         1.0e12)) {
                std::cerr << "error edge bound del_hat_global in ze" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_bound.GetID()][3 * dof] << ' '
                          << rhs_global_prev[edge_bound.GetID()][3 * dof] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_bound.GetID()][3 * dof + 1],
                                         rhs_global_prev[edge_bound.GetID()][3 * dof + 1],
                                         1.0e12)) {
                std::cerr << "error edge bound del_hat_global in qx" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_bound.GetID()][3 * dof + 1] << ' '
                          << rhs_global_prev[edge_bound.GetID()][3 * dof + 1] << std::endl;

                any_error = true;
            }

            if (!Utilities::almost_equal(delta_global_diff_est[edge_bound.GetID()][3 * dof + 2],
                                         rhs_global_prev[edge_bound.GetID()][3 * dof + 2],
                                         1.0e12)) {
                std::cerr << "error edge bound del_hat_global in qy" << std::endl;
                std::cout << std::setprecision(15) << delta_global_diff_est[edge_bound.GetID()][3 * dof + 2] << ' '
                          << rhs_global_prev[edge_bound.GetID()][3 * dof + 2] << std::endl;

                any_error = true;
            }
        }
    });

    if (any_error) {
        return 1;
    }

    return 0;
}