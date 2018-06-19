#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/lcos.hpp>

#include "../../../general_definitions.hpp"
#include "../swe_definitions.hpp"

#include "../function_files/swe_initial_condition_functions.hpp"
#include "../function_files/swe_source_functions.hpp"
#include "../function_files/swe_true_solution_functions.hpp"

#include "../swe_problem.hpp"
#include "../kernels_preprocessor/swe_kernels_preprocessor.hpp"
#include "../kernels_processor/swe_kernels_processor.hpp"
#include "../kernels_postprocessor/swe_kernels_postprocessor.hpp"

#include "../../../simulation/rkdg_simulation/rkdg_simulation_hpx.hpp"
#include "../../../simulation/stepper/rk_stepper.hpp"

using hpx_simulation_unit_swe           = RKDG::HPXSimulationUnit<SWE::Problem>;
using hpx_simulation_unit_swe_component = hpx::components::simple_component<RKDG::HPXSimulationUnit<SWE::Problem>>;
HPX_REGISTER_COMPONENT(hpx_simulation_unit_swe_component, hpx_simulation_unit_swe);

using hpx_simulation_swe           = RKDG::HPXSimulation<SWE::Problem>;
using hpx_simulation_swe_component = hpx::components::simple_component<RKDG::HPXSimulation<SWE::Problem>>;
HPX_REGISTER_COMPONENT(hpx_simulation_swe_component, hpx_simulation_swe);

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

    std::vector<RKDG::HPXSimulationClient<SWE::Problem>> simulation_clients;
    simulation_clients.reserve(localities.size());

    auto t1 = std::chrono::high_resolution_clock::now();
    for (hpx::naming::id_type const& locality : localities) {
        hpx::future<hpx::id_type> simulation_id =
            hpx::new_<hpx::components::simple_component<RKDG::HPXSimulation<SWE::Problem>>>(locality, input_string);

        simulation_clients.emplace_back(std::move(simulation_id));
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

    return hpx::finalize();
}