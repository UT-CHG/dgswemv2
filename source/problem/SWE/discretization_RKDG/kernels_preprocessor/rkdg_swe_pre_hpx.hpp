#ifndef RKDG_SWE_PRE_HPX_HPP
#define RKDG_SWE_PRE_HPX_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace RKDG {
template <typename HPXSimUnitType>
auto Problem::preprocessor_hpx(HPXSimUnitType* sim_unit) {
    initialize_data_parallel_pre_send(sim_unit->discretization.mesh, sim_unit->problem_input, CommTypes::baryctr_coord);

    hpx::future<void> receive_future =
        sim_unit->communicator.ReceiveAll(CommTypes::baryctr_coord, sim_unit->stepper.GetTimestamp());

    sim_unit->communicator.SendAll(CommTypes::baryctr_coord, sim_unit->stepper.GetTimestamp());

    return receive_future.then([sim_unit](auto&&) {
        initialize_data_parallel_post_receive(sim_unit->discretization.mesh, CommTypes::baryctr_coord);

        sim_unit->discretization.mesh.CallForEachElement(
            [sim_unit](auto& elt) { elt.data.resize(sim_unit->stepper.GetNumStages() + 1); });
    });
}
}
}

#endif