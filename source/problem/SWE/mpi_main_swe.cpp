#include <mpi.h>

#include "general_definitions.hpp"

#include "swe_problem.hpp"
#include "swe_kernels_preprocessor.hpp"
#include "swe_kernels_processor.hpp"
#include "swe_kernels_postprocessor.hpp"

#include "simulation/mpi_simulation.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/DG_HYPER_SWE input_file\n";
        return 1;
    } else {
        auto t1 = std::chrono::high_resolution_clock::now();
        MPI_Init(NULL, NULL);

        std::string input_string = std::string(argv[1]);

        MPISimulation<SWE::Problem> simulation(input_string);

        simulation.Run();

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Finalize();
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        return 0;
    }
}
