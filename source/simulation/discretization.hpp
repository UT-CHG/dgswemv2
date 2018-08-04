#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"

template <typename ProblemType>
struct DGDiscretization {
    typename ProblemType::ProblemMeshType mesh;

    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input, Writer<ProblemType>& writer) {
        std::tuple<> empty_comm;

        initialize_mesh<ProblemType>(this->mesh, input, empty_comm, writer);

        ProblemType::initialize_data_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);
    }

    template <typename CommunicatorType>
    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input,
                    CommunicatorType& communicator,
                    Writer<ProblemType>& writer) {
        initialize_mesh<ProblemType>(this->mesh, input, communicator, writer);

        ProblemType::initialize_data_parallel_pre_send_kernel(
            this->mesh, input.mesh_input.mesh_data, input.problem_input);
    }
};

template <typename ProblemType>
struct HDGDiscretization {
    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemMeshSkeletonType mesh_skeleton;

    /* have to template on global data structure class */
    SparseMatrix<double> delta_local_inv;
    SparseMatrix<double> delta_hat_local;
    DynVector<double> rhs_local;

    SparseMatrix<double> delta_global;
    SparseMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;

    SparseMatrix<double> w1_w1;
    SparseMatrix<double> w1_w2;
    SparseMatrix<double> w1_w1_hat;
    DynVector<double> w1_rhs;

    SparseMatrix<double> w2_w1;
    SparseMatrix<double> w2_w2_inv;
    SparseMatrix<double> w2_w1_hat;

    SparseMatrix<double> w1_hat_w1;
    SparseMatrix<double> w1_hat_w2;
    SparseMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;
    /* have to template on global data structure class */

    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input, Writer<ProblemType>& writer) {
        std::tuple<> empty_comm;

        initialize_mesh<ProblemType>(this->mesh, input, empty_comm, writer);
        initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, writer);

        ProblemType::initialize_data_kernel(this->mesh, input.mesh_input.mesh_data, input.problem_input);

        ProblemType::initialize_global_problem(this);
    }

    template <typename CommunicatorType>
    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input,
                    CommunicatorType& communicator,
                    Writer<ProblemType>& writer) {
        initialize_mesh<ProblemType>(this->mesh, input, communicator, writer);
        initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, writer);

        ProblemType::initialize_data_parallel_pre_send_kernel(
            this->mesh, input.mesh_input.mesh_data, input.problem_input);

        ProblemType::initialize_global_problem(this);
    }
};

#endif