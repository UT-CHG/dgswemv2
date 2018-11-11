#include "general_definitions.hpp"

#include "simulation/serial/simulation_base.hpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/dgswemv2-serial input_file\n";
        return 1;
    } else {
        std::string input_string = std::string(argv[1]);

        std::unique_ptr<Serial::SimulationBase> simulation = Serial::SimulationFactory::Create(input_string);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulation->Run();
        auto t2 = std::chrono::high_resolution_clock::now();

        std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
                  << "\n";

        simulation->Finalize();

        return 0;
    }
}
