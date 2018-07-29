#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/lcos.hpp>

#include "general_definitions.hpp"
#include "problem/GN/gn_definitions.hpp"

#include "problem/GN/problem_function_files/gn_initial_condition_functions.hpp"
#include "problem/GN/problem_function_files/gn_source_functions.hpp"
#include "problem/GN/problem_function_files/gn_true_solution_functions.hpp"

#include "problem/GN/discretization_EHDG/ehdg_gn_problem.hpp"
#include "problem/GN/discretization_EHDG/kernels_preprocessor/ehdg_gn_kernels_preprocessor.hpp"
#include "problem/GN/discretization_EHDG/kernels_processor/ehdg_gn_proc_hpx_stage.hpp"
#include "problem/GN/discretization_EHDG/kernels_postprocessor/ehdg_gn_kernels_postprocessor.hpp"

#include "simulation/hpx/simulation_hpx.hpp"
#include "simulation/stepper/rk_stepper.hpp"

REGISTER_HPX_COMPONENTS(GN::EHDG::Problem);

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/EHDG_GN_HPX input_file\n";
        return 1;
    } else {
        return hpx::init(argc, argv);
    }
}

int hpx_main(int argc, char* argv[]) {
    std::string input_string = std::string(argv[1]);

    const std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

    std::vector<HPXSimulationClient<GN::EHDG::Problem>> simulation_clients;
    simulation_clients.reserve(localities.size());

    auto t1 = std::chrono::high_resolution_clock::now();
    for (hpx::naming::id_type const& locality : localities) {
        hpx::future<hpx::id_type> simulation_id =
            hpx::new_<hpx::components::simple_component<HPXSimulation<GN::EHDG::Problem>>>(locality, input_string);

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
