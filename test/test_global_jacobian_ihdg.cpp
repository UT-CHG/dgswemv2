#include "general_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_IHDG/ihdg_swe_problem.hpp"
#include "problem/SWE/discretization_IHDG/kernels_preprocessor/ihdg_swe_kernels_preprocessor.hpp"
#include "problem/SWE/discretization_IHDG/kernels_processor/ihdg_swe_kernels_processor.hpp"
#include "problem/SWE/discretization_IHDG/kernels_postprocessor/ihdg_swe_kernels_postprocessor.hpp"

#include "simulation/simulation_IHDG/serial/ihdg_simulation.hpp"
#include "simulation/stepper/rk_stepper.hpp"

#include "utilities/almost_equal.hpp"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    std::string input_string = "../../test/files_for_testing/jacobians/input.15";

    InputParameters<typename SWE::IHDG::Problem::ProblemInputType> input(input_string);

    typename SWE::IHDG::Problem::ProblemMeshType mesh;
    typename SWE::IHDG::Problem::ProblemMeshSkeletonType mesh_skeleton;

    RKStepper stepper;
    Writer<SWE::IHDG::Problem> writer;

    input.read_mesh();  // read mesh meta data
    input.read_bcis();  // read bc data

    mesh = typename SWE::IHDG::Problem::ProblemMeshType(input.polynomial_order);

    stepper = RKStepper(input.stepper_input);
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

    SWE::IHDG::Problem::initialize_data_kernel(mesh, input.mesh_input.mesh_data, input.problem_input);

    mesh.CallForEachElement([](auto& elt) { elt.data.resize(2); });

    DynMatrix<double> delta_local;
    DynMatrix<double> delta_local_inv;
    DynMatrix<double> delta_hat_local;
    DynVector<double> rhs_local;

    DynMatrix<double> delta_global;
    DynMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;

    uint dof_local = mesh.GetNumberElements() * 9;  // hardcode for p=1
    uint dof_global =
        mesh_skeleton.GetNumberEdgeInterfaces() * 6 + mesh_skeleton.GetNumberEdgeBoundaries() * 6;  // hardcode for p=1

    delta_local.resize(dof_local, dof_local);
    delta_local_inv.resize(dof_local, dof_local);
    delta_hat_local.resize(dof_local, dof_global);
    rhs_local.resize(dof_local);

    delta_global.resize(dof_global, dof_local);
    delta_hat_global.resize(dof_global, dof_global);
    rhs_global.resize(dof_global);

    delta_local     = 0.0;
    delta_local_inv = 0.0;
    delta_hat_local = 0.0;
    rhs_local       = 0.0;

    delta_global     = 0.0;
    delta_hat_global = 0.0;
    rhs_global       = 0.0;

    mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        // Set IDs for global matrix construction
        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            elt.data.boundary[bound_id].elt_ID = elt.GetID();
        }

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
    });

    mesh_skeleton.CallForEachEdgeInterface([](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // Set IDs for global matrix construction
        boundary_in.edg_ID = edge_int.GetID();
        boundary_ex.edg_ID = edge_int.GetID();

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

    mesh_skeleton.CallForEachEdgeBoundary([](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // Set IDs for global matrix construction
        boundary.edg_ID = edge_bound.GetID();

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

        for (uint dof = 0; dof < elt.data.get_ndof(); dof++) {
            state.q[dof][SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q[dof][SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            state.q[dof][SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); dof++) {
            edge_state.q_hat[dof][SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat[dof][SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat[dof][SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); dof++) {
            edge_state.q_hat[dof][SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat[dof][SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            edge_state.q_hat[dof][SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
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

    // assemble global
    mesh.CallForEachElement([&](auto& elt) {
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        submatrix(delta_local, elt_ID * 9, elt_ID * 9, 9, 9)     = internal.delta_local;
        submatrix(delta_local_inv, elt_ID * 9, elt_ID * 9, 9, 9) = inverse(internal.delta_local);
        subvector(rhs_local, elt_ID * 9, 9)                      = internal.rhs_local;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            uint edg_ID = elt.data.boundary[bound_id].edg_ID;

            submatrix(delta_hat_local, elt_ID * 9, edg_ID * 6, 9, 6) = elt.data.boundary[bound_id].delta_hat_local;
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_in_ID * 9, 6, 9) = boundary_in.delta_global;
        submatrix(delta_global, edg_ID * 6, elt_ex_ID * 9, 6, 9) = boundary_ex.delta_global;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_ID * 9, 6, 9) = boundary.delta_global;
    });

    // delat q and q_hats to be assigned randomly
    DynVector<double> delta_q(dof_local);
    DynVector<double> delta_q_hat(dof_global);

    // randomly assign delta_q
    for (uint dof = 0; dof < dof_local; dof++) {
        delta_q[dof] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
    }

    delta_q *= 1.0e-8;  // make it small

    // randomly assign delta_q_hat
    for (uint dof = 0; dof < dof_global; dof++) {
        delta_q_hat[dof] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
    }

    delta_q_hat *= 1.0e-8;  // make it small

    // containers to store data
    DynVector<double> rhs_local_prev(dof_local);
    DynVector<double> delta_hat_local_diff_est(dof_local);
    DynVector<double> delta_local_diff_est(dof_local);
    // containers to store data

    // store for future comparison
    rhs_local_prev = rhs_local;

    // estimate difference with jacobian
    delta_hat_local_diff_est = delta_hat_local * delta_q_hat;

    // containers to store data
    DynVector<double> rhs_global_prev(dof_global);
    DynVector<double> delta_hat_global_diff_est(dof_global);
    DynVector<double> delta_global_diff_est(dof_global);
    // containers to store data

    // store for future comparison
    rhs_global_prev = rhs_global;

    // estimate difference with jacobian
    delta_hat_global_diff_est = delta_hat_global * delta_q_hat;

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        uint edg_ID = edge_int.GetID();

        auto rhs_global_ref = subvector(delta_q_hat, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global_ref[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global_ref[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global_ref[3 * dof + 2];
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        uint edg_ID = edge_bound.GetID();

        auto rhs_global_ref = subvector(delta_q_hat, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global_ref[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global_ref[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global_ref[3 * dof + 2];
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

    // assemble global
    mesh.CallForEachElement([&](auto& elt) {
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        submatrix(delta_local, elt_ID * 9, elt_ID * 9, 9, 9)     = internal.delta_local;
        submatrix(delta_local_inv, elt_ID * 9, elt_ID * 9, 9, 9) = inverse(internal.delta_local);
        subvector(rhs_local, elt_ID * 9, 9)                      = internal.rhs_local;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            uint edg_ID = elt.data.boundary[bound_id].edg_ID;

            submatrix(delta_hat_local, elt_ID * 9, edg_ID * 6, 9, 6) = elt.data.boundary[bound_id].delta_hat_local;
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_in_ID * 9, 6, 9) = boundary_in.delta_global;
        submatrix(delta_global, edg_ID * 6, elt_ex_ID * 9, 6, 9) = boundary_ex.delta_global;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_ID * 9, 6, 9) = boundary.delta_global;
    });

    // get true diff
    rhs_local_prev -= rhs_local;
    rhs_global_prev -= rhs_global;

    for (uint dof = 0; dof < dof_local; dof++) {
        if (!Utilities::almost_equal(delta_hat_local_diff_est[dof], rhs_local_prev[dof], 1.0e12)) {
            std::cerr << "error del hat local" << std::endl;
            std::cout << std::setprecision(15) << delta_hat_local_diff_est[dof] << ' ' << rhs_local_prev[dof]
                      << std::endl;
        }
    }

    for (uint dof = 0; dof < dof_global; dof++) {
        if (!Utilities::almost_equal(delta_hat_global_diff_est[dof], rhs_global_prev[dof], 1.0e12)) {
            std::cerr << "error del hat global" << std::endl;
            std::cout << std::setprecision(15) << delta_hat_global_diff_est[dof] << ' ' << rhs_global_prev[dof]
                      << std::endl;
        }
    }

    rhs_local_prev       = rhs_local;
    delta_local_diff_est = delta_local * delta_q;

    rhs_global_prev       = rhs_global;
    delta_global_diff_est = delta_global * delta_q;

    mesh.CallForEachElement([&](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        uint elt_ID = elt.GetID();

        auto rhs_local_ref = subvector(delta_q, elt_ID * 9, 9);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q[dof][SWE::Variables::ze] += rhs_local_ref[3 * dof];
            state.q[dof][SWE::Variables::qx] += rhs_local_ref[3 * dof + 1];
            state.q[dof][SWE::Variables::qy] += rhs_local_ref[3 * dof + 2];
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

    // assemble global
    mesh.CallForEachElement([&](auto& elt) {
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        submatrix(delta_local, elt_ID * 9, elt_ID * 9, 9, 9)     = internal.delta_local;
        submatrix(delta_local_inv, elt_ID * 9, elt_ID * 9, 9, 9) = inverse(internal.delta_local);
        subvector(rhs_local, elt_ID * 9, 9)                      = internal.rhs_local;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            uint edg_ID = elt.data.boundary[bound_id].edg_ID;

            submatrix(delta_hat_local, elt_ID * 9, edg_ID * 6, 9, 6) = elt.data.boundary[bound_id].delta_hat_local;
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_in_ID * 9, 6, 9) = boundary_in.delta_global;
        submatrix(delta_global, edg_ID * 6, elt_ex_ID * 9, 6, 9) = boundary_ex.delta_global;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_ID * 9, 6, 9) = boundary.delta_global;
    });

    // get true diff
    rhs_local_prev -= rhs_local;
    rhs_global_prev -= rhs_global;

    for (uint dof = 0; dof < dof_local; dof++) {
        if (!Utilities::almost_equal(delta_local_diff_est[dof], rhs_local_prev[dof], 1.0e12)) {
            std::cerr << "error del local" << std::endl;
            std::cout << std::setprecision(15) << delta_local_diff_est[dof] << ' ' << rhs_local_prev[dof] << std::endl;
        }
    }

    for (uint dof = 0; dof < dof_global; dof++) {
        if (!Utilities::almost_equal(delta_global_diff_est[dof], rhs_global_prev[dof], 1.0e12)) {
            std::cerr << "error del global" << std::endl;
            std::cout << std::setprecision(15) << delta_global_diff_est[dof] << ' ' << rhs_global_prev[dof]
                      << std::endl;
        }
    }

    // randomly assign delta_q
    for (uint dof = 0; dof < dof_local; dof++) {
        delta_q[dof] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
    }

    delta_q *= 1.0e-8;  // make it small

    // randomly assign delta_q_hat
    for (uint dof = 0; dof < dof_global; dof++) {
        delta_q_hat[dof] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
    }

    delta_q_hat *= 1.0e-8;  // make it small

    DynMatrix<double> global;
    DynVector<double> rhs;
    DynVector<double> delta;

    DynVector<double> diff_est;
    DynVector<double> rhs_prev;

    global.resize(dof_local + dof_global, dof_local + dof_global);
    rhs.resize(dof_local + dof_global);
    delta.resize(dof_local + dof_global);

    submatrix(global, 0, 0, dof_local, dof_local)                   = delta_local;
    submatrix(global, 0, dof_local, dof_local, dof_global)          = delta_hat_local;
    submatrix(global, dof_local, 0, dof_global, dof_local)          = delta_global;
    submatrix(global, dof_local, dof_local, dof_global, dof_global) = delta_hat_global;

    subvector(rhs, 0, dof_local)          = rhs_local;
    subvector(rhs, dof_local, dof_global) = rhs_global;

    subvector(delta, 0, dof_local)          = delta_q;
    subvector(delta, dof_local, dof_global) = delta_q_hat;

    rhs_prev = rhs;
    diff_est = global * delta;

    mesh.CallForEachElement([&](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        uint elt_ID = elt.GetID();

        auto rhs_local_ref = subvector(delta_q, elt_ID * 9, 9);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q[dof][SWE::Variables::ze] += rhs_local_ref[3 * dof];
            state.q[dof][SWE::Variables::qx] += rhs_local_ref[3 * dof + 1];
            state.q[dof][SWE::Variables::qy] += rhs_local_ref[3 * dof + 2];
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        uint edg_ID = edge_int.GetID();

        auto rhs_global_ref = subvector(delta_q_hat, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global_ref[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global_ref[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global_ref[3 * dof + 2];
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        uint edg_ID = edge_bound.GetID();

        auto rhs_global_ref = subvector(delta_q_hat, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global_ref[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global_ref[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global_ref[3 * dof + 2];
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

    // assemble global
    mesh.CallForEachElement([&](auto& elt) {
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        submatrix(delta_local, elt_ID * 9, elt_ID * 9, 9, 9)     = internal.delta_local;
        submatrix(delta_local_inv, elt_ID * 9, elt_ID * 9, 9, 9) = inverse(internal.delta_local);
        subvector(rhs_local, elt_ID * 9, 9)                      = internal.rhs_local;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            uint edg_ID = elt.data.boundary[bound_id].edg_ID;

            submatrix(delta_hat_local, elt_ID * 9, edg_ID * 6, 9, 6) = elt.data.boundary[bound_id].delta_hat_local;
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_in_ID * 9, 6, 9) = boundary_in.delta_global;
        submatrix(delta_global, edg_ID * 6, elt_ex_ID * 9, 6, 9) = boundary_ex.delta_global;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_ID * 9, 6, 9) = boundary.delta_global;
    });

    subvector(rhs, 0, dof_local)          = rhs_local;
    subvector(rhs, dof_local, dof_global) = rhs_global;

    // get true diff
    rhs_prev -= rhs;

    for (uint dof = 0; dof < dof_local + dof_global; dof++) {
        if (!Utilities::almost_equal(diff_est[dof], rhs_prev[dof], 1.0e12)) {
            std::cerr << "error del compound" << std::endl;
            std::cout << std::setprecision(15) << diff_est[dof] << ' ' << rhs_prev[dof] << std::endl;
        }
    }

    // randomly assign delta_q_hat
    for (uint dof = 0; dof < dof_global; dof++) {
        delta_q_hat[dof] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
    }

    delta_q_hat *= 1.0e-8;  // make it small

    delta_q = delta_local_inv * (rhs_local - delta_local_inv * delta_hat_local * delta_q_hat);

    delta_q *= 1.0e-8; // delta_q is too large cause we have delta_local_inv * rhs_local term in there

    rhs_global_prev       = rhs_global;
    delta_global_diff_est = delta_hat_global * delta_q_hat + delta_global * delta_q;

    mesh.CallForEachElement([&](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        uint elt_ID = elt.GetID();

        auto rhs_local_ref = subvector(delta_q, elt_ID * 9, 9);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q[dof][SWE::Variables::ze] += rhs_local_ref[3 * dof];
            state.q[dof][SWE::Variables::qx] += rhs_local_ref[3 * dof + 1];
            state.q[dof][SWE::Variables::qy] += rhs_local_ref[3 * dof + 2];
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        uint edg_ID = edge_int.GetID();

        auto rhs_global_ref = subvector(delta_q_hat, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global_ref[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global_ref[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global_ref[3 * dof + 2];
        }
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        uint edg_ID = edge_bound.GetID();

        auto rhs_global_ref = subvector(delta_q_hat, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global_ref[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global_ref[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global_ref[3 * dof + 2];
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

    // assemble global
    mesh.CallForEachElement([&](auto& elt) {
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        submatrix(delta_local_inv, elt_ID * 9, elt_ID * 9, 9, 9) = inverse(internal.delta_local);
        subvector(rhs_local, elt_ID * 9, 9)                      = internal.rhs_local;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            uint edg_ID = elt.data.boundary[bound_id].edg_ID;

            submatrix(delta_hat_local, elt_ID * 9, edg_ID * 6, 9, 6) = elt.data.boundary[bound_id].delta_hat_local;
        }
    });

    mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_in_ID * 9, 6, 9) = boundary_in.delta_global;
        submatrix(delta_global, edg_ID * 6, elt_ex_ID * 9, 6, 9) = boundary_ex.delta_global;
    });

    mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        submatrix(delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        subvector(rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        submatrix(delta_global, edg_ID * 6, elt_ID * 9, 6, 9) = boundary.delta_global;
    });

    rhs_global_prev -= rhs_global;

    for (uint dof = 0; dof < dof_global; dof++) {
        if (!Utilities::almost_equal(delta_global_diff_est[dof], rhs_global_prev[dof], 1.0e12)) {
            std::cerr << "error back substitute del global" << std::endl;
            std::cout << std::setprecision(15) << delta_global_diff_est[dof] << ' ' << rhs_global_prev[dof]
                      << std::endl;
        }
    }
}
