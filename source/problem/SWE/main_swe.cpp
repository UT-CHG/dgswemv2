#include "../../general_definitions.hpp"

#include "../../preprocessor/input_parameters.hpp"
#include "../../stepper.hpp"
#include "../../initialize_mesh.hpp"
#include "../../run_simulation.hpp"

#include "swe_problem.hpp"
#include "swe_kernels.hpp"

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>

void local_main(std::string);
HPX_PLAIN_ACTION(local_main, local_main_act);

void solve_mesh(std::string);
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
    
    std::vector<hpx::future<void> > futures;
    futures.reserve(localities.size());

    for (hpx::naming::id_type const& node : localities) {
        futures.push_back(hpx::async<local_main_act>(node, std::string(argv[1])));
    }

    hpx::wait_all(futures);

    return hpx::finalize();  // Handles HPX shutdown
}

void local_main(std::string input_string) {
    const uint n_threads = hpx::get_os_thread_count();
    const hpx::naming::id_type here = hpx::find_here();

    std::vector<hpx::future<void> > futures;
    futures.reserve(n_threads);

    for(uint thread = 0; thread < n_threads; thread++){
        futures.push_back(hpx::async<solve_mesh_act>(here, input_string));
    }

    hpx::wait_all(futures);
}

void solve_mesh(std::string input_string) {
    hpx::cout << input_string << '\n';
    /*
    try {
        const InputParameters input(input_string.c_str());

        printf("Starting program %s with p=%d for %s mesh\n\n",
               argv[0],
               input.polynomial_order,
               input.mesh_file_name.c_str());

        auto mesh = initialize_mesh<SWE::Problem>(input.polynomial_order, input.mesh_data);

        SWE::Problem::initialize_data_kernel(*mesh, input.mesh_data);
        
        Stepper stepper(input.rk.nstages, input.rk.order, input.dt);

        auto t1 = std::chrono::high_resolution_clock::now();
        run_simulation<SWE::Problem>(input.T_end, stepper, *mesh);
        auto t2 = std::chrono::high_resolution_clock::now();
        
        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Exception caught\n";
        std::cerr << "  " << e.what() << std::endl;
    }*/
}
