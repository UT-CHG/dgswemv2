#ifndef IHDG_SWE_PRE_OMPI_HPP
#define IHDG_SWE_PRE_OMPI_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace IHDG {
template <template <typename> class OMPISimUnitType, typename ProblemType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                                typename ProblemType::ProblemGlobalDataType& global_data,
                                const ProblemStepperType& stepper,
                                const uint begin_sim_id,
                                const uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::baryctr_coord, 0);

        initialize_data_parallel_pre_send(
            sim_units[su_id]->discretization.mesh, sim_units[su_id]->problem_input, CommTypes::baryctr_coord);

        sim_units[su_id]->communicator.SendAll(CommTypes::baryctr_coord, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::baryctr_coord, 0);

        initialize_data_parallel_post_receive(sim_units[su_id]->discretization.mesh, CommTypes::baryctr_coord);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::baryctr_coord, 0);
    }

#pragma omp barrier
#pragma omp master
    {
        std::vector<uint> global_dof_offsets;

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            uint global_dof_offset = 0;

            Problem::initialize_global_problem_parallel_pre_send(
                sim_units[su_id]->discretization, stepper, global_dof_offset);

            global_dof_offsets.push_back(global_dof_offset);
        }

        uint total_global_dof_offset = 0;

        if (!global_dof_offsets.empty()) {
            total_global_dof_offset = std::accumulate(global_dof_offsets.begin(), global_dof_offsets.end(), 0);

            // exclusive scan
            std::rotate(global_dof_offsets.begin(), global_dof_offsets.end() - 1, global_dof_offsets.end());

            global_dof_offsets.front() = 0;

            for (uint su_id = 1; su_id < sim_units.size(); ++su_id) {
                global_dof_offsets[su_id] += global_dof_offsets[su_id - 1];
            }
        }

        int n_localities;
        int locality_id;

        MPI_Comm_size(MPI_COMM_WORLD, &n_localities);
        MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

        std::vector<uint> total_global_dof_offsets;

        if (locality_id == 0) {
            total_global_dof_offsets.resize(n_localities);
        }

        MPI_Gather(&total_global_dof_offset,
                   1,
                   MPI_UNSIGNED,
                   &total_global_dof_offsets.front(),
                   1,
                   MPI_UNSIGNED,
                   0,
                   MPI_COMM_WORLD);

        uint n_global_dofs;

        if (locality_id == 0) {
            n_global_dofs = std::accumulate(total_global_dof_offsets.begin(), total_global_dof_offsets.end(), 0);

            // exclusive scan
            std::rotate(
                total_global_dof_offsets.begin(), total_global_dof_offsets.end() - 1, total_global_dof_offsets.end());

            total_global_dof_offsets.front() = 0;

            for (int locality_id = 1; locality_id < n_localities; ++locality_id) {
                total_global_dof_offsets[locality_id] += total_global_dof_offsets[locality_id - 1];
            }
        }

        MPI_Bcast(&n_global_dofs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        MPI_Scatter(&total_global_dof_offsets.front(),
                    1,
                    MPI_UNSIGNED,
                    &total_global_dof_offset,
                    1,
                    MPI_UNSIGNED,
                    0,
                    MPI_COMM_WORLD);

        MatCreate(MPI_COMM_WORLD, &(global_data.delta_hat_global));
        MatSetSizes(global_data.delta_hat_global, PETSC_DECIDE, PETSC_DECIDE, n_global_dofs, n_global_dofs);
        MatSetUp(global_data.delta_hat_global);

        VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, n_global_dofs, &(global_data.rhs_global));

        KSPCreate(MPI_COMM_WORLD, &(global_data.ksp));
        KSPSetOperators(global_data.ksp, global_data.delta_hat_global, global_data.delta_hat_global);

        KSPGetPC(global_data.ksp, &(global_data.pc));
        PCSetType(global_data.pc, PCLU);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->communicator.ReceiveAll(CommTypes::init_global_prob, 0);

            global_dof_offsets[su_id] += total_global_dof_offset;

            Problem::initialize_global_problem_parallel_finalize_pre_send(sim_units[su_id]->discretization,
                                                                          global_dof_offsets[su_id]);

            sim_units[su_id]->communicator.SendAll(CommTypes::init_global_prob, 0);
        }

        std::vector<uint> global_dof_indx;
        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->communicator.WaitAllReceives(CommTypes::init_global_prob, 0);

            Problem::initialize_global_problem_parallel_post_receive(sim_units[su_id]->discretization, global_dof_indx);
        }

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->communicator.WaitAllSends(CommTypes::init_global_prob, 0);
        }

        VecCreateSeq(MPI_COMM_SELF, global_dof_indx.size(), &(global_data.sol));

        ISCreateGeneral(MPI_COMM_SELF,
                        global_dof_indx.size(),
                        (int*)&global_dof_indx.front(),
                        PETSC_COPY_VALUES,
                        &(global_data.from));
        ISCreateStride(MPI_COMM_SELF, global_dof_indx.size(), 0, 1, &(global_data.to));

        VecScatterCreate(
            global_data.rhs_global, global_data.from, global_data.sol, global_data.to, &(global_data.scatter));
    }
#pragma omp barrier

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { elt.data.resize(stepper.GetNumStages() + 1); });
    }
}
}
}

#endif
