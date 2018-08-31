#include "general_definitions.hpp"
#include "problem/Green-Naghdi/gn_definitions.hpp"

#include "problem/Green-Naghdi/problem_function_files/gn_initial_condition_functions.hpp"
#include "problem/Green-Naghdi/problem_function_files/gn_source_functions.hpp"
#include "problem/Green-Naghdi/problem_function_files/gn_true_solution_functions.hpp"

#include "problem/Green-Naghdi/discretization_IHDG/ihdg_gn_problem.hpp"
#include "problem/Green-Naghdi/discretization_IHDG/kernels_preprocessor/ihdg_gn_kernels_preprocessor.hpp"
#include "problem/Green-Naghdi/discretization_IHDG/kernels_postprocessor/ihdg_gn_kernels_postprocessor.hpp"

#include "problem/Green-Naghdi/discretization_IHDG/kernels_preprocessor/ihdg_gn_pre_serial.hpp"
#include "problem/Green-Naghdi/discretization_IHDG/kernels_processor/ihdg_gn_proc_serial_stage.hpp"

#include "simulation/serial/simulation.hpp"
#include "simulation/stepper/rk_stepper.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/ihdg_GN_SERIAL input_file\n";
        return 1;
    } else {
        std::string input_string = std::string(argv[1]);

        Simulation<GN::IHDG::Problem> simulation(input_string);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulation.Run();
        auto t2 = std::chrono::high_resolution_clock::now();

        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        return 0;
    }
}
