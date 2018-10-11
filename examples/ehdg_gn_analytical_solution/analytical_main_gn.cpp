#include "general_definitions.hpp"

// SWE part of the problem
#include "problem/SWE/swe_definitions.hpp"

#include "analytical_gn_initial_condition_functions.hpp"
#include "analytical_gn_source_functions.hpp"
#include "analytical_gn_true_solution_functions.hpp"

#include "problem/SWE/discretization_EHDG/ehdg_swe_problem.hpp"
#include "problem/SWE/discretization_EHDG/kernels_preprocessor/ehdg_swe_kernels_preprocessor.hpp"
#include "problem/SWE/discretization_EHDG/kernels_postprocessor/ehdg_swe_kernels_postprocessor.hpp"

#include "problem/SWE/discretization_EHDG/kernels_preprocessor/ehdg_swe_pre_serial.hpp"
#include "problem/SWE/discretization_EHDG/kernels_processor/ehdg_swe_proc_serial_stage.hpp"

// GN part of the problem
#include "problem/Green-Naghdi/gn_definitions.hpp"

#include "problem/Green-Naghdi/discretization_EHDG/ehdg_gn_problem.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/kernels_preprocessor/ehdg_gn_kernels_preprocessor.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/kernels_postprocessor/ehdg_gn_kernels_postprocessor.hpp"

#include "problem/Green-Naghdi/discretization_EHDG/kernels_preprocessor/ehdg_gn_pre_serial.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/kernels_processor/ehdg_gn_proc_serial_stage.hpp"

#include "simulation/serial/simulation.hpp"
#include "simulation/stepper/rk_stepper.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/bin input_file\n";
        return 1;
    } else {
        std::string input_string = std::string(argv[1]);

        Simulation<GN::EHDG::Problem> simulation(input_string);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulation.Run();
        auto t2 = std::chrono::high_resolution_clock::now();

        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        simulation.ComputeL2Residual();

        return 0;
    }
}
