#ifndef IHDG_SWE_PRE_OMPI_HPP
#define IHDG_SWE_PRE_OMPI_HPP

#include "ihdg_swe_pre_init_data.hpp"

namespace SWE {
namespace IHDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    std::vector<uint> local_dof_offsets;
    std::vector<uint> global_dof_offsets;

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        Problem::initialize_data_parallel(sim_units[su_id]->discretization.mesh, sim_units[su_id]->problem_input);

        uint local_dof_offset  = 0;
        uint global_dof_offset = 0;

        Problem::initialize_global_problem_parallel_pre_send(
            sim_units[su_id]->discretization, sim_units[su_id]->communicator, local_dof_offset, global_dof_offset);

        local_dof_offsets.push_back(local_dof_offset);
        global_dof_offsets.push_back(global_dof_offset);
    }

    uint total_local_dof_offset  = std::accumulate(local_dof_offsets.begin(), local_dof_offsets.end(), 0);
    uint total_global_dof_offset = std::accumulate(global_dof_offsets.begin(), global_dof_offsets.end(), 0);

    // exclusive scan
    std::rotate(local_dof_offsets.begin(), local_dof_offsets.end() - 1, local_dof_offsets.end());
    std::rotate(global_dof_offsets.begin(), global_dof_offsets.end() - 1, global_dof_offsets.end());

    local_dof_offsets.front()  = 0;
    global_dof_offsets.front() = 0;

    for (uint su_id = 1; su_id < sim_units.size(); ++su_id) {
        local_dof_offsets[su_id] += local_dof_offsets[su_id - 1];
        global_dof_offsets[su_id] += global_dof_offsets[su_id - 1];
    }

    int n_localities;
    int locality_id;

    MPI_Comm_size(MPI_COMM_WORLD, &n_localities);
    MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

    std::vector<uint> total_local_dof_offsets;
    std::vector<uint> total_global_dof_offsets;

    if (locality_id == 0) {
        total_local_dof_offsets.resize(n_localities);
        total_global_dof_offsets.resize(n_localities);
    }

    MPI_Gather(
        &total_local_dof_offset, 1, MPI_UNSIGNED, &total_local_dof_offsets.front(), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    MPI_Gather(&total_global_dof_offset,
               1,
               MPI_UNSIGNED,
               &total_global_dof_offsets.front(),
               1,
               MPI_UNSIGNED,
               0,
               MPI_COMM_WORLD);

    if (locality_id == 0) {
        // exclusive scan
        std::rotate(total_local_dof_offsets.begin(), total_local_dof_offsets.end() - 1, total_local_dof_offsets.end());
        std::rotate(
            total_global_dof_offsets.begin(), total_global_dof_offsets.end() - 1, total_global_dof_offsets.end());

        total_local_dof_offsets.front()  = 0;
        total_global_dof_offsets.front() = 0;

        for (int locality_id = 1; locality_id < n_localities; ++locality_id) {
            total_local_dof_offsets[locality_id] += total_local_dof_offsets[locality_id - 1];
            total_global_dof_offsets[locality_id] += total_global_dof_offsets[locality_id - 1];
        }
    }

    MPI_Scatter(
        &total_local_dof_offsets.front(), 1, MPI_UNSIGNED, &total_local_dof_offset, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    MPI_Scatter(&total_global_dof_offsets.front(),
                1,
                MPI_UNSIGNED,
                &total_global_dof_offset,
                1,
                MPI_UNSIGNED,
                0,
                MPI_COMM_WORLD);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        local_dof_offsets[su_id] += total_local_dof_offset;
        global_dof_offsets[su_id] += total_global_dof_offset;
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        Problem::initialize_global_problem_parallel_finalize_pre_send(
            sim_units[su_id]->discretization, local_dof_offsets[su_id], global_dof_offsets[su_id]);
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(SWE::CommTypes::preprocessor,
                                                  sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.SendAll(SWE::CommTypes::preprocessor, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(SWE::CommTypes::preprocessor,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        Problem::initialize_global_problem_parallel_post_receive(sim_units[su_id]->discretization,
                                                                 sim_units[su_id]->communicator);
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(SWE::CommTypes::preprocessor,
                                                    sim_units[su_id]->stepper.GetTimestamp());
    }
}
}
}

#endif