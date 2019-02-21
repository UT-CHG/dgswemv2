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

    uint dof_local  = 0;
    uint dof_global = 0;

    discretization.mesh.CallForEachElement([&](auto& elt) {
        elt.data.resize(2);

        // randomly assign q
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage + 1];

        dof_local += elt.data.get_ndof() * SWE::n_variables;

        // randomly assign q
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                state.q(var, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        dof_global += edge_int.edge_data.get_ndof() * SWE::n_variables;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                edge_state.q_hat(var, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        dof_global += edge_bound.edge_data.get_ndof() * SWE::n_variables;

        // randomly assign q_hat
        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                edge_state.q_hat(var, dof) = -1.0 + 2.0 * ((double)rand() / (RAND_MAX));
            }
        }
    });

    // do one pass to compute all jacobians and rhs
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
    // do one pass to compute all jacobians and rhs

    SparseMatrix<double>& delta_hat_global = global_data.delta_hat_global;
    DynVector<double>& rhs_global          = global_data.rhs_global;

    SparseMatrixMeta<double> sparse_delta_hat_global;

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        internal.delta_local_inv = inverse(internal.delta_local);
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global, &sparse_delta_hat_global](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        edge_internal.delta_hat_global -=
            boundary_in.delta_global * internal_in.delta_local_inv * boundary_in.delta_hat_local +
            boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_ex.delta_hat_local;

        edge_internal.rhs_global -= boundary_in.delta_global * internal_in.delta_local_inv * internal_in.rhs_local +
                                    boundary_ex.delta_global * internal_ex.delta_local_inv * internal_ex.rhs_local;

        subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size()) = edge_internal.rhs_global;

        for (uint i = 0; i < global_dof_indx.size(); ++i) {
            for (uint j = 0; j < global_dof_indx.size(); ++j) {
                sparse_delta_hat_global.add_triplet(
                    global_dof_indx[i], global_dof_indx[j], edge_internal.delta_hat_global(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_in)
                continue;

            auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary_in.delta_global * internal_in.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_ex)
                continue;

            auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global, &sparse_delta_hat_global](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        edge_internal.delta_hat_global -= boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

        edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

        subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size()) = edge_internal.rhs_global;

        for (uint i = 0; i < global_dof_indx.size(); ++i) {
            for (uint j = 0; j < global_dof_indx.size(); ++j) {
                sparse_delta_hat_global.add_triplet(
                    global_dof_indx[i], global_dof_indx[j], edge_internal.delta_hat_global(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
            if (bound_id == edge_bound.boundary.bound_id)
                continue;

            auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }
    });

    sparse_delta_hat_global.get_sparse_matrix(delta_hat_global);

    // randomly assign delta_q_hat
    DynVector<double> delta_q_hat(dof_global);

    for (uint dof = 0; dof < dof_global; ++dof) {
        delta_q_hat[dof] = 1.0e-8 * (-1.0 + 2.0 * ((double)rand() / (RAND_MAX)));
    }

    // containers to store data
    DynVector<double> delta_hat_diff(dof_global);
    DynVector<double> delta_hat_diff_est(dof_global);

    // store for future comparison
    delta_hat_diff = rhs_global;

    // estimate difference with jacobian
    delta_hat_diff_est = delta_hat_global * delta_q_hat;

    discretization.mesh_skeleton.CallForEachEdgeInterface([&delta_q_hat](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        auto dq_hat = subvector(delta_q_hat, (uint)global_dof_indx[0], (uint)global_dof_indx.size());

        internal_in.rhs_local -= boundary_in.delta_hat_local * dq_hat;
        internal_ex.rhs_local -= boundary_ex.delta_hat_local * dq_hat;

        edge_state.q_hat += reshape<double, SWE::n_variables, SO::ColumnMajor>(dq_hat, edge_int.edge_data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&delta_q_hat](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        auto dq_hat = subvector(delta_q_hat, (uint)global_dof_indx[0], (uint)global_dof_indx.size());

        internal.rhs_local -= boundary.delta_hat_local * dq_hat;

        edge_state.q_hat += reshape<double, SWE::n_variables, SO::ColumnMajor>(dq_hat, edge_bound.edge_data.get_ndof());
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        auto& internal = elt.data.internal;

        internal.rhs_local = internal.delta_local_inv * internal.rhs_local;

        internal.rhs_local *= 1.0e-8;  // delta_q is too large cause we have delta_local_inv * rhs_local term in there

        state.q += reshape<double, SWE::n_variables, SO::ColumnMajor>(internal.rhs_local, elt.data.get_ndof());
    });

    // do one pass to compute all jacobians and rhs
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
    // do one pass to compute all jacobians and rhs

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        internal.delta_local_inv = inverse(internal.delta_local);
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global, &sparse_delta_hat_global](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        edge_internal.delta_hat_global -=
            boundary_in.delta_global * internal_in.delta_local_inv * boundary_in.delta_hat_local +
            boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_ex.delta_hat_local;

        edge_internal.rhs_global -= boundary_in.delta_global * internal_in.delta_local_inv * internal_in.rhs_local +
                                    boundary_ex.delta_global * internal_ex.delta_local_inv * internal_ex.rhs_local;

        subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size()) = edge_internal.rhs_global;

        for (uint i = 0; i < global_dof_indx.size(); ++i) {
            for (uint j = 0; j < global_dof_indx.size(); ++j) {
                sparse_delta_hat_global.add_triplet(
                    global_dof_indx[i], global_dof_indx[j], edge_internal.delta_hat_global(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_in)
                continue;

            auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary_in.delta_global * internal_in.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_ex)
                continue;

            auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global, &sparse_delta_hat_global](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        edge_internal.delta_hat_global -= boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

        edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

        subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size()) = edge_internal.rhs_global;

        for (uint i = 0; i < global_dof_indx.size(); ++i) {
            for (uint j = 0; j < global_dof_indx.size(); ++j) {
                sparse_delta_hat_global.add_triplet(
                    global_dof_indx[i], global_dof_indx[j], edge_internal.delta_hat_global(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
            if (bound_id == edge_bound.boundary.bound_id)
                continue;

            auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }
    });

    // get true diff
    delta_hat_diff -= rhs_global;

    for (uint dof = 0; dof < dof_global; ++dof) {
        if (!Utilities::almost_equal(delta_hat_diff_est[dof], delta_hat_diff[dof], 1.0e12)) {
            std::cerr << "error del hat global" << std::endl;
            std::cout << std::setprecision(15) << delta_hat_diff_est[dof] << ' ' << delta_hat_diff[dof] << std::endl;

            any_error = true;
        }
    }

    //if (any_error) {
    //    return 1;
    //}

    return 0;
}
