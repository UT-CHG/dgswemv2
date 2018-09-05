#ifndef EHDG_SWE_PRE_HPX_HPP
#define EHDG_SWE_PRE_HPX_HPP

#include "ehdg_swe_pre_init_data.hpp"

namespace SWE {
namespace EHDG {
template <typename HPXSimUnitType>
decltype(auto) Problem::preprocessor_hpx(HPXSimUnitType* sim_unit) {
    Problem::initialize_data_parallel(sim_unit->discretization.mesh, sim_unit->problem_input);

    Problem::initialize_global_problem(sim_unit->discretization);

    return hpx::make_ready_future();
}
}
}

#endif