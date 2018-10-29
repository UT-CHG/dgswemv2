#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/lcos.hpp>

#include "general_definitions.hpp"
#include "problem/SWE/swe_definitions.hpp"

#include "manufactured_swe_initial_condition_functions.hpp"
#include "manufactured_swe_source_functions.hpp"
#include "manufactured_swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"
#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_kernels_preprocessor.hpp"

#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_pre_hpx.hpp"
#include "problem/SWE/discretization_RKDG/kernels_processor/rkdg_swe_proc_hpx_stage.hpp"

#include "simulation/hpx/simulation_hpx.hpp"
#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"

REGISTER_HPX_COMPONENTS(SWE::RKDG::Problem);

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/DG_HYPER_SWE input_file\n";
        return 1;
    } else {
        return hpx::init(argc, argv);
    }
}

int hpx_main(int argc, char* argv[]) {
    std::string input_string = std::string(argv[1]);

    const std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

    std::vector<HPXSimulationClient<SWE::RKDG::Problem>> simulation_clients;
    simulation_clients.reserve(localities.size());

    auto t1 = std::chrono::high_resolution_clock::now();
    for (hpx::naming::id_type const& locality : localities) {
        simulation_clients.emplace_back(hpx::new_<HPXSimulation<SWE::RKDG::Problem>>(locality, input_string));
    }

    std::vector<hpx::future<void>> run_futures;
    run_futures.reserve(simulation_clients.size());

    for (auto& sim_client : simulation_clients) {
        run_futures.push_back(sim_client.Run());
    }

    hpx::wait_all(run_futures);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << std::endl;

    hpx::future<double> globalResidualL2 = ComputeL2Residual(simulation_clients);
    std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(globalResidualL2.get()) << std::endl;

    return hpx::finalize();
}
