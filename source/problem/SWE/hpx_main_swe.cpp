#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>

#include "../../general_definitions.hpp"

#include "../../preprocessor/input_parameters.hpp"
#include "../../stepper.hpp"
#include "../../preprocessor/initialize_mesh.hpp"

#include "swe_problem.hpp"
#include "swe_kernels.hpp"

#include "../../simulation/hpx_simulation.hpp"

// using hpx_mesh_swe_component = hpx::components::simple_component<hpx_mesh<SWE::Problem>>;
// using hpx_mesh_swe = hpx_mesh<SWE::Problem>;
// HPX_REGISTER_COMPONENT(hpx_mesh_swe_component, hpx_mesh_swe);

void local_main(std::string);
HPX_PLAIN_ACTION(local_main, local_main_act);

hpx::future<void> solve_mesh(std::string, uint);
HPX_PLAIN_ACTION(solve_mesh, solve_mesh_act);

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
    const std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

    std::vector<hpx::future<void>> futures;
    futures.reserve(localities.size());

    auto t1 = std::chrono::high_resolution_clock::now();
    for (hpx::naming::id_type const& node : localities) {
        futures.push_back(hpx::async<local_main_act>(node, std::string(argv[1])));
    }

    hpx::wait_all(futures);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << "\n";

    return hpx::finalize();  // Handles HPX shutdown
}

void local_main(std::string input_string) {
    const uint n_threads = hpx::get_os_thread_count();
    const hpx::naming::id_type here = hpx::find_here();

    std::vector<hpx::future<void>> futures;
    futures.reserve(n_threads);

    for (uint thread = 0; thread < n_threads; thread++) {
        futures.push_back(hpx::async<solve_mesh_act>(here, input_string, thread));
    }

    hpx::wait_all(futures);
}

hpx::future<void> solve_mesh(std::string input_string, uint thread) {
    try {
        HPXSimulation<SWE::Problem> simulation(input_string, hpx::get_locality_id(), thread);

        return simulation.Run(43200.0);
    }
    catch (const std::exception& e) {
        std::cerr << "Exception caught\n";
        std::cerr << "  " << e.what() << std::endl;
    }
}
