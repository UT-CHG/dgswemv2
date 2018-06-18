#include "../../../general_definitions.hpp"
#include "../swe_definitions.hpp"

#include "../function_files/swe_initial_condition_functions.hpp"
#include "../function_files/swe_source_functions.hpp"
#include "../function_files/swe_true_solution_functions.hpp"

#include "../swe_problem.hpp"
#include "../kernels_preprocessor/swe_kernels_preprocessor.hpp"
#include "../kernels_processor/swe_kernels_processor.hpp"
#include "../kernels_postprocessor/swe_kernels_postprocessor.hpp"

#include "../../../simulation/rkdg_simulation/rkdg_simulation.hpp"
#include "../../../simulation/rkdg_simulation/rkdg_stepper.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/DG_HYPER_SWE input_file\n";
        return 1;
    } else {
        std::string input_string = std::string(argv[1]);

        RKDGSimulation<SWE::Problem> simulation(input_string);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulation.Run();
        auto t2 = std::chrono::high_resolution_clock::now();

        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        return 0;
    }
}
