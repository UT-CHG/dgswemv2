#include <mpi.h>
#include <omp.h>

#include "general_definitions.hpp"
#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"
#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_kernels_preprocessor.hpp"

#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_pre_ompi.hpp"
#include "problem/SWE/discretization_RKDG/kernels_processor/rkdg_swe_proc_ompi_step.hpp"

#include "simulation/ompi/simulation_ompi.hpp"
#include "simulation/stepper/rk_stepper.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/RKDG_SWE_OMPI input_file\n";
        return 1;
    } else {
        auto t1 = std::chrono::high_resolution_clock::now();
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

        if (provided != MPI_THREAD_MULTIPLE) {
            if (omp_get_max_threads() > 1) {
                std::cerr << "dgswemv2 with ompi parallelization was submitted with more than \n"
                          << "1 thread per MPI rank and MPI_THREAD_MULTIPLE is not provided!\n"
                          << "Please either find an MPI implementation that supports MPI_THREAD_MULTIPLE,\n"
                          << " or resubmit the job with one thread per MPI rank and set the\n"
                          << " environment variable OMP_NUM_THREADS=1\n";
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
                return 1;
            }
        }

        std::string input_string = std::string(argv[1]);

        OMPISimulation<SWE::RKDG::Problem> simulation(input_string);

        simulation.Run();

        MPI_Barrier(MPI_COMM_WORLD);
        auto t2 = std::chrono::high_resolution_clock::now();

        int locality_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

        if (locality_id == 0) {
            std::cout << "Time Elapsed (in us): "
                      << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
        }

        MPI_Finalize();

        return 0;
    }
}
