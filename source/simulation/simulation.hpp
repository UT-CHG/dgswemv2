#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"

template <typename ProblemType>
class Simulation {
  protected:
    InputParameters<typename ProblemType::InputType> input;
    Stepper stepper;
    typename ProblemType::ProblemMeshType mesh;
    std::ofstream log_file;

  public:
    Simulation() : input(), stepper(input.rk.nstages, input.rk.order, input.dt),
            mesh(input.polynomial_order) {}

    Simulation(std::string& input_string) :
            input(input_string),
            stepper(input.rk.nstages, input.rk.order, input.dt),
            mesh(input.polynomial_order) {
        input.ReadMesh();
        mesh.SetMeshName(input.mesh_data.mesh_name);
        log_file = std::ofstream("output/" + input.mesh_data.mesh_name + "_log",
                std::ofstream::out);
        log_file << "Starting simulation with p=" << input.polynomial_order << " for "
                 << mesh.GetMeshName() << " mesh" << std::endl << std::endl;
        std::tuple<> empty_comm;
        initialize_mesh<ProblemType>(this->mesh, input.mesh_data, empty_comm, input.problem_input);
    }

    ~Simulation() {
      log_file.close();
    }
};

#endif
