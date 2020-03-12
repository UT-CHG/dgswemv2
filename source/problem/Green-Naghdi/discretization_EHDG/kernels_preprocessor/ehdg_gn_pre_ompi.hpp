#ifndef EHDG_GN_PRE_OMPI_HPP
#define EHDG_GN_PRE_OMPI_HPP

#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_init_data.hpp"
#include "ehdg_gn_pre_dbath_ompi.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                ProblemGlobalDataType& global_data,
                                const ProblemStepperType& stepper,
                                const uint begin_sim_id,
                                const uint end_sim_id) {
    SWE_SIM::Problem::preprocessor_ompi(sim_units, global_data, stepper.GetFirstStepper(), begin_sim_id, end_sim_id);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        GN::initialize_data_parallel_pre_send(sim_units[su_id]->discretization.mesh, SWE_SIM::CommTypes::baryctr_coord);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(SWE_SIM::CommTypes::baryctr_coord, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.SendAll(SWE_SIM::CommTypes::baryctr_coord, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(SWE_SIM::CommTypes::baryctr_coord, 0);

        GN::initialize_data_parallel_post_receive(sim_units[su_id]->discretization.mesh,
                                                  SWE_SIM::CommTypes::baryctr_coord);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(SWE_SIM::CommTypes::baryctr_coord, 0);
    }

#pragma omp barrier
#pragma omp master
    {
        uint dc_global_dof_offset = 0;
        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            Problem::initialize_global_dc_problem_parallel_pre_send(sim_units[su_id]->discretization,
                                                                    dc_global_dof_offset);
        }

        int n_localities;
        int locality_id;
        MPI_Comm_size(MPI_COMM_WORLD, &n_localities);
        MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

        std::vector<uint> dc_global_dof_offsets;
        if (locality_id == 0) {
            dc_global_dof_offsets.resize(n_localities);
        }
        MPI_Gather(
            &dc_global_dof_offset, 1, MPI_UNSIGNED, &dc_global_dof_offsets.front(), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        uint n_dc_global_dofs = 0;
        if (locality_id == 0) {
            n_dc_global_dofs = std::accumulate(dc_global_dof_offsets.begin(), dc_global_dof_offsets.end(), 0);
        }
        MPI_Bcast(&n_dc_global_dofs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        const uint n_global_dofs = (sim_units.front()->discretization.mesh.GetP() + 1) * GN::n_dimensions;
        MatCreateBAIJ(MPI_COMM_WORLD,
                      n_global_dofs,
                      n_global_dofs * dc_global_dof_offset,
                      n_global_dofs * dc_global_dof_offset,
                      n_global_dofs * n_dc_global_dofs,
                      n_global_dofs * n_dc_global_dofs,
                      5,
                      NULL,
                      4,
                      NULL,
                      &(global_data.w1_hat_w1_hat));

        VecCreateMPI(MPI_COMM_WORLD,
                     n_global_dofs * dc_global_dof_offset,
                     n_global_dofs * n_dc_global_dofs,
                     &(global_data.w1_hat_rhs));
        VecSetBlockSize(global_data.w1_hat_rhs, n_global_dofs);

        KSPCreate(MPI_COMM_WORLD, &(global_data.dc_ksp));
        KSPSetOperators(global_data.dc_ksp, global_data.w1_hat_w1_hat, global_data.w1_hat_w1_hat);
        // KSPGetPC(global_data.dc_ksp, &(global_data.dc_pc));
        // PCSetType(global_data.dc_pc, PCLU);

        PetscLogStageRegister("Construct", &global_data.con_stage);
        PetscLogStageRegister("Solve", &global_data.sol_stage);
        PetscLogStageRegister("Propagate", &global_data.prop_stage);
        PetscLogStageRegister("SWE", &global_data.swe_stage);
        PetscLogStageRegister("Derivatives", &global_data.d_stage);

        if (locality_id == 0) {
            std::rotate(dc_global_dof_offsets.begin(), dc_global_dof_offsets.end() - 1, dc_global_dof_offsets.end());
            dc_global_dof_offsets.front() = 0;
            for (int locality_id = 1; locality_id < n_localities; ++locality_id) {
                dc_global_dof_offsets[locality_id] += dc_global_dof_offsets[locality_id - 1];
            }
        }
        MPI_Scatter(
            &dc_global_dof_offsets.front(), 1, MPI_UNSIGNED, &dc_global_dof_offset, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->communicator.ReceiveAll(CommTypes::dc_global_dof_indx, 0);
        }

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            Problem::initialize_global_dc_problem_parallel_finalize_pre_send(sim_units[su_id]->discretization,
                                                                             dc_global_dof_offset);
            sim_units[su_id]->communicator.SendAll(CommTypes::dc_global_dof_indx, 0);
        }

        std::vector<uint> dc_global_dof_indx;
        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dc_global_dof_indx, 0);
            Problem::initialize_global_dc_problem_parallel_post_receive(sim_units[su_id]->discretization,
                                                                        dc_global_dof_indx);
        }

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->communicator.WaitAllSends(CommTypes::dc_global_dof_indx, 0);
        }

        VecCreateSeq(MPI_COMM_SELF, n_global_dofs * dc_global_dof_indx.size(), &(global_data.dc_sol));

        ISCreateBlock(MPI_COMM_SELF,
                      n_global_dofs,
                      dc_global_dof_indx.size(),
                      (int*)&dc_global_dof_indx.front(),
                      PETSC_COPY_VALUES,
                      &(global_data.dc_from));
        ISCreateStride(MPI_COMM_SELF, n_global_dofs * dc_global_dof_indx.size(), 0, 1, &(global_data.dc_to));

        VecScatterCreate(global_data.w1_hat_rhs,
                         global_data.dc_from,
                         global_data.dc_sol,
                         global_data.dc_to,
                         &(global_data.dc_scatter));

#if defined(B_RECON_AVG) || defined(D_RECON_AVG)
        std::set<uint> nodeIDs;
        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&nodeIDs](auto& elt) { nodeIDs.insert(elt.GetNodeID().begin(), elt.GetNodeID().end()); });
        }
        global_data.local_nodeIDs = std::vector<int>(nodeIDs.begin(), nodeIDs.end());
        std::map<uint, uint> nodeIDmap;
        uint it = 0;
        for (auto nodeID : nodeIDs) {
            nodeIDmap[nodeID] = it++;
        }

        uint local_max_nodeID  = *std::max_element(nodeIDs.begin(), nodeIDs.end());
        uint global_max_nodeID = 0;
        MPI_Allreduce(&local_max_nodeID, &global_max_nodeID, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

        VecCreateMPI(MPI_COMM_WORLD,
                     (global_max_nodeID + 1) * GN::n_dddbath_terms,
                     PETSC_DECIDE,
                     &(global_data.global_derivatives_at_node));
        VecSetBlockSize(global_data.global_derivatives_at_node, GN::n_dddbath_terms);
        VecCreateSeq(MPI_COMM_SELF,
                     global_data.local_nodeIDs.size() * GN::n_dddbath_terms,
                     &(global_data.local_derivatives_at_node));
        global_data.derivatives_at_node = DynVector<double>(global_data.local_nodeIDs.size() * GN::n_dddbath_terms);

        ISCreateBlock(MPI_COMM_SELF,
                      GN::n_dddbath_terms,
                      global_data.local_nodeIDs.size(),
                      (int*)&global_data.local_nodeIDs.front(),
                      PETSC_COPY_VALUES,
                      &(global_data.derivatives_from));
        ISCreateStride(
            MPI_COMM_SELF, global_data.local_nodeIDs.size() * GN::n_dddbath_terms, 0, 1, &(global_data.derivatives_to));
        VecScatterCreate(global_data.global_derivatives_at_node,
                         global_data.derivatives_from,
                         global_data.local_derivatives_at_node,
                         global_data.derivatives_to,
                         &(global_data.derivatives_scatter));

        VecCreateMPI(MPI_COMM_WORLD, global_max_nodeID + 1, PETSC_DECIDE, &(global_data.global_node_mult));
        VecCreateSeq(MPI_COMM_SELF, global_data.local_nodeIDs.size(), &(global_data.local_node_mult));

        ISCreateGeneral(MPI_COMM_SELF,
                        global_data.local_nodeIDs.size(),
                        (int*)&global_data.local_nodeIDs.front(),
                        PETSC_COPY_VALUES,
                        &(global_data.node_mult_from));
        ISCreateStride(MPI_COMM_SELF, global_data.local_nodeIDs.size(), 0, 1, &(global_data.node_mult_to));
        VecScatterCreate(global_data.global_node_mult,
                         global_data.node_mult_from,
                         global_data.local_node_mult,
                         global_data.node_mult_to,
                         &(global_data.node_mult_scatter));

        VecSet(global_data.global_node_mult, 0.0);
        DynVector<double> node_mult(global_data.local_nodeIDs.size());
        set_constant(node_mult, 0.0);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([&node_mult, &nodeIDmap](auto& elt) {
                auto& derivative = elt.data.derivative;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    derivative.local_nodeID[node] = nodeIDmap[elt.GetNodeID()[node]];
                    node_mult[derivative.local_nodeID[node]] += 1.0;
                }
            });
        }

        VecSetValues(global_data.global_node_mult,
                     global_data.local_nodeIDs.size(),
                     &global_data.local_nodeIDs.front(),
                     node_mult.data(),
                     ADD_VALUES);

        VecAssemblyBegin(global_data.global_node_mult);
        VecAssemblyEnd(global_data.global_node_mult);

        VecScatterBegin(global_data.node_mult_scatter,
                        global_data.global_node_mult,
                        global_data.local_node_mult,
                        INSERT_VALUES,
                        SCATTER_FORWARD);
        VecScatterEnd(global_data.node_mult_scatter,
                      global_data.global_node_mult,
                      global_data.local_node_mult,
                      INSERT_VALUES,
                      SCATTER_FORWARD);

        double* m_ptr;
        VecGetArray(global_data.local_node_mult, &m_ptr);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([m_ptr](auto& elt) {
                auto& derivative = elt.data.derivative;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    derivative.node_mult[node] = m_ptr[derivative.local_nodeID[node]];
                }
            });
        }

        VecRestoreArray(global_data.local_node_mult, &m_ptr);
#endif
    }
#pragma omp barrier

    Problem::compute_bathymetry_derivatives_ompi(sim_units, global_data, begin_sim_id, end_sim_id);

    uint n_stages = stepper.GetFirstStepper().GetNumStages() > stepper.GetSecondStepper().GetNumStages()
                        ? stepper.GetFirstStepper().GetNumStages()
                        : stepper.GetSecondStepper().GetNumStages();

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [n_stages](auto& elt) { elt.data.resize(n_stages + 1); });
    }
}
}
}

#endif
