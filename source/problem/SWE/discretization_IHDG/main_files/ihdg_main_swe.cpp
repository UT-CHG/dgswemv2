#include "general_definitions.hpp"
#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_IHDG/ihdg_swe_problem.hpp"
#include "problem/SWE/discretization_IHDG/kernels_preprocessor/ihdg_swe_kernels_preprocessor.hpp"

#include "problem/SWE/discretization_IHDG/kernels_preprocessor/ihdg_swe_pre_serial.hpp"
#include "problem/SWE/discretization_IHDG/kernels_processor/ihdg_swe_proc_serial_step.hpp"

#include "simulation/serial/simulation.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/IHDG_SWE_SERIAL input_file\n";
        return 1;
    } else {
        std::string input_string = std::string(argv[1]);

        Simulation<SWE::IHDG::Problem> simulation(input_string);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulation.Run();
        auto t2 = std::chrono::high_resolution_clock::now();

        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        return 0;
    }
}
