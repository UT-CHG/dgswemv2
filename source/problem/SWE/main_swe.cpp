#include "../../general_definitions.hpp"

#include "../../preprocessor/input_parameters.hpp"

#include "../../stepper.hpp"

#include "../../initialize_mesh.hpp"
#include "../../run_simulation.hpp"

#include "swe_problem.hpp"
#include "swe_kernels.hpp"

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/DG_HYPER_SWE input_file\n";
        return 1;
    }

    try {
        const InputParameters input(argv[1]);

        printf("Starting program %s with p=%d for %s mesh\n\n",
               argv[1],
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

        delete mesh;

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Exception caught\n";
        std::cerr << "  " << e.what() << std::endl;

        return 1;
    }
}