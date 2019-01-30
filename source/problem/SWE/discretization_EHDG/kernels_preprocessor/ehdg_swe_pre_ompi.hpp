#ifndef EHDG_SWE_PRE_OMPI_HPP
#define EHDG_SWE_PRE_OMPI_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                uint begin_sim_id,
                                uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        initialize_data_parallel_pre_send(
            sim_units[su_id]->discretization.mesh, sim_units[su_id]->problem_input, CommTypes::baryctr_coord);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::baryctr_coord, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.SendAll(CommTypes::baryctr_coord, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::baryctr_coord,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        initialize_data_parallel_post_receive(sim_units[su_id]->discretization.mesh, CommTypes::baryctr_coord);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::baryctr_coord, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        Problem::initialize_global_problem_parallel_pre_send(sim_units[su_id]->discretization);
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::init_global_prob,
                                                  sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.SendAll(CommTypes::init_global_prob, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::init_global_prob,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        Problem::initialize_global_problem_parallel_post_receive(sim_units[su_id]->discretization);
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::init_global_prob,
                                                    sim_units[su_id]->stepper.GetTimestamp());
    }
}
}
}

#endif