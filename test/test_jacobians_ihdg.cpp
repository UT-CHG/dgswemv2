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

    SparseMatrix<double> delta_local_inv;
    SparseMatrix<double> delta_hat_local;
    DynVector<double> rhs_local;

    SparseMatrix<double> delta_global;
    SparseMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;

    uint dof_local = mesh.GetNumberElements() * 9;  // hardcode for p=1
    uint dof_global =
        mesh_skeleton.GetNumberEdgeInterfaces() * 6 + mesh_skeleton.GetNumberEdgeBoundaries() * 6;  // hardcode for p=1

    delta_local_inv.resize(dof_local, dof_local);
    delta_hat_local.resize(dof_local, dof_global);
    rhs_local.resize(dof_local);

    delta_global.resize(dof_global, dof_local);
    delta_hat_global.resize(dof_global, dof_global);
    rhs_global.resize(dof_global);

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
        auto& state      = elt.data.state[stage];

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

    mesh.CallForEachElement([&](auto& elt) { SWE::IHDG::Problem::local_volume_kernel(stepper, elt); });

    mesh.CallForEachInterface([&](auto& intface) { SWE::IHDG::Problem::local_interface_kernel(stepper, intface); });

    mesh.CallForEachBoundary([&](auto& bound) { SWE::IHDG::Problem::local_boundary_kernel(stepper, bound); });

    mesh_skeleton.CallForEachEdgeInterface(
        [&](auto& edge_int) { SWE::IHDG::Problem::local_edge_interface_kernel(stepper, edge_int); });

    mesh_skeleton.CallForEachEdgeBoundary(
        [&](auto& edge_bound) { SWE::IHDG::Problem::local_edge_boundary_kernel(stepper, edge_bound); });

    /*mesh.CallForEachElement([&](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];

        // randomly assign delta_q and increment
        DynVector<double> delta_q(SWE::n_variables * elt.data.get_ndof());

        for (uint dof = 0; dof < elt.data.get_ndof(); dof++) {
            delta_q[3 * dof + SWE::Variables::ze] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q[3 * dof + SWE::Variables::qx] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            delta_q[3 * dof + SWE::Variables::qy] = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
        }

        delta_q *= 1.0e-8;  // make it small

        for (uint dof = 0; dof < elt.data.get_ndof(); dof++) {
            state.q[dof][SWE::Variables::ze] += delta_q[3 * dof];
            state.q[dof][SWE::Variables::qx] += delta_q[3 * dof + 1];
            state.q[dof][SWE::Variables::qy] += delta_q[3 * dof + 2];
        }
        // randomly assign delta_q and increment

        // difference estimate
        DynVector<double> diff_est(SWE::n_variables * elt.data.get_ndof());

        diff_est = elt.data.internal.delta_local * delta_q;
        // difference estimate

        // subtract from rhs for future comparison
        elt.data.internal.rhs_local =- diff_est;
    });*/
}