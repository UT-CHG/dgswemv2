#include <mpi.h>
#include <omp.h>
#ifdef HAS_PETSC
#include <petscksp.h>
#endif

#include "general_definitions.hpp"

#include "simulation/ompi/simulation_ompi_base.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/dgswemv2-ompi input_file\n";
        return 1;
    } else {
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

#ifdef HAS_PETSC
        PetscInitialize(&argc, &argv, (char*)0, NULL);
#endif

        std::string input_string = std::string(argv[1]);

        std::unique_ptr<OMPISimulationBase> simulation = OMPISimulationFactory::Create(input_string);

	simulation->Initialize();
	MPI_Barrier(MPI_COMM_WORLD);
	auto t1 = std::chrono::high_resolution_clock::now();

        simulation->Run();

        MPI_Barrier(MPI_COMM_WORLD);
        auto t2 = std::chrono::high_resolution_clock::now();

        int locality_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

        if (locality_id == 0) {
            std::cout << "Time Elapsed (in us): "
                      << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
        }

        simulation->Finalize();

#ifdef HAS_PETSC
        PetscFinalize();
#endif

        MPI_Finalize();

        return 0;
    }
}
