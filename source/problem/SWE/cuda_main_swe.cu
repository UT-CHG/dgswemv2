#include "../../general_definitions.hpp"

#include "simulation/stepper.hpp"

#include "swe_problem.hpp"
#include "swe_kernels_preprocessor.hpp"
#include "swe_cuda_kernels_processor.hpp"
#include "swe_kernels_postprocessor.hpp"

#include "simulation/cuda_simulation.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/DG_HYPER_SWE input_file\n";
        return 1;
    } else {
        std::string input_string = std::string(argv[1]);

        CUDASimulation<SWE::CUDAProblem> simulation(input_string);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulation.Run();
        auto t2 = std::chrono::high_resolution_clock::now();

        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        return 0;
    }
}
