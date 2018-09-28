#ifndef RKDG_SWE_PRE_OMPI_HPP
#define RKDG_SWE_PRE_OMPI_HPP

#include "rkdg_swe_pre_init_data.hpp"

namespace SWE {
namespace RKDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        Problem::initialize_data_parallel_pre_send(sim_units[su_id]->discretization.mesh,
                                                   sim_units[su_id]->problem_input);
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::baryctr_coord,
                                                  sim_units[su_id]->stepper.GetTimestamp());
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.SendAll(CommTypes::baryctr_coord,
                                               sim_units[su_id]->stepper.GetTimestamp());
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::baryctr_coord,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        Problem::initialize_data_parallel_post_receive(sim_units[su_id]->discretization.mesh);
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::baryctr_coord,
                                                    sim_units[su_id]->stepper.GetTimestamp());
    }
}
}
}

#endif