#ifndef EHDG_SWE_PRE_OMPI_HPP
#define EHDG_SWE_PRE_OMPI_HPP

namespace SWE {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                uint begin_sim_id,
                                uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        Problem::initialize_data_parallel(sim_units[su_id]->discretization.mesh, sim_units[su_id]->problem_input);

        Problem::initialize_global_problem(sim_units[su_id]->discretization);
    }
}
}
}

#endif