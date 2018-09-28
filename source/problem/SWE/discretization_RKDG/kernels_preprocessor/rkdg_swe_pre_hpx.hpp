#ifndef RKDG_SWE_PRE_HPX_HPP
#define RKDG_SWE_PRE_HPX_HPP

#include "rkdg_swe_pre_init_data.hpp"

namespace SWE {
namespace RKDG {
template <typename HPXSimUnitType>
decltype(auto) Problem::preprocessor_hpx(HPXSimUnitType* sim_unit) {
    Problem::initialize_data_parallel_pre_send(sim_unit->discretization.mesh, sim_unit->problem_input);

    hpx::future<void> receive_future =
        sim_unit->communicator.ReceiveAll(SWE::RKDG::CommTypes::baryctr_coord, sim_unit->stepper.GetTimestamp());

    sim_unit->communicator.SendAll(SWE::RKDG::CommTypes::baryctr_coord, sim_unit->stepper.GetTimestamp());

    return receive_future.then(
        [sim_unit](auto&&) { Problem::initialize_data_parallel_post_receive(sim_unit->discretization.mesh); });
}
}
}

#endif