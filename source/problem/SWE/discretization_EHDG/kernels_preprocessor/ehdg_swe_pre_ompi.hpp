#ifndef EHDG_SWE_PRE_OMPI_HPP
#define EHDG_SWE_PRE_OMPI_HPP

namespace SWE {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        Problem::initialize_data_parallel(sim_units[su_id]->discretization.mesh, sim_units[su_id]->problem_input);

        Problem::initialize_global_problem(sim_units[su_id]->discretization);
    }
}
}
}

#endif