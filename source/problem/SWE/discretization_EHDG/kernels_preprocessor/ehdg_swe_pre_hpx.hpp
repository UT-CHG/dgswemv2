#ifndef EHDG_SWE_PRE_HPX_HPP
#define EHDG_SWE_PRE_HPX_HPP

#include "ehdg_swe_pre_init_data.hpp"

namespace SWE {
namespace EHDG {
template <typename HPXSimUnitType>
decltype(auto) Problem::hpx_preprocessor_kernel(HPXSimUnitType* sim_unit) {
    Problem::initialize_data_parallel_pre_send_kernel(sim_unit->discretization.mesh, sim_unit->problem_input);

    hpx::future<void> receive_future =
        sim_unit->communicator.ReceiveAll(SWE::CommTypes::preprocessor, sim_unit->stepper.GetTimestamp());

    sim_unit->communicator.SendAll(SWE::CommTypes::preprocessor, sim_unit->stepper.GetTimestamp());

    return receive_future.then(
        [sim_unit](auto&&) { Problem::initialize_data_parallel_post_receive_kernel(sim_unit->discretization.mesh); });
}
}
}

#endif