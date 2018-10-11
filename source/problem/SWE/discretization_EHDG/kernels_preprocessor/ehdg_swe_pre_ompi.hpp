#ifndef EHDG_SWE_PRE_OMPI_HPP
#define EHDG_SWE_PRE_OMPI_HPP

namespace SWE {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

    n_threads = (uint)omp_get_num_threads();
    thread_id = (uint)omp_get_thread_num();

    sim_per_thread = (sim_units.size() + n_threads - 1) / n_threads;

    begin_sim_id = sim_per_thread * thread_id;
    end_sim_id   = std::min(sim_per_thread * (thread_id + 1), (uint)sim_units.size());

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        Problem::initialize_data_parallel(sim_units[su_id]->discretization.mesh, sim_units[su_id]->problem_input);

        Problem::initialize_global_problem(sim_units[su_id]->discretization);
    }
}
}
}

#endif