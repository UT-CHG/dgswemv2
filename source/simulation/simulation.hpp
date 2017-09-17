#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"

template <typename ProblemType>
class Simulation {
  private:
    InputParameters input;

    Stepper stepper;
    typename ProblemType::ProblemMeshType mesh;

    std::string log_file_name;

  public:
    Simulation() : input(), stepper(input.rk.nstages, input.rk.order, input.dt), mesh(input.polynomial_order) {}
    Simulation(std::string input_string)
        : input(input_string), stepper(input.rk.nstages, input.rk.order, input.dt), mesh(input.polynomial_order) {
        input.ReadMesh();

        mesh.SetMeshName(input.mesh_data.mesh_name);

        this->log_file_name = "output/" + input.mesh_data.mesh_name + "_log";

        std::ofstream log_file(this->log_file_name, std::ofstream::out);

        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << mesh.GetMeshName() << " mesh"
                 << std::endl << std::endl;

        std::tuple<> empty_comm;

        initialize_mesh<ProblemType>(this->mesh, input.mesh_data, empty_comm);
    }

    void Run();
};

template <typename ProblemType>
void Simulation<ProblemType>::Run() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    auto scrutinize_solution_kernel = [this](auto& elt) {
        ProblemType::scrutinize_solution_kernel(this->stepper, elt);
    };

    uint nsteps = (uint)std::ceil(this->input.T_end / this->stepper.get_dt());
    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    ProblemType::write_modal_data_kernel(this->stepper, this->mesh);

    for (uint step = 1; step <= nsteps; ++step) {
        for (uint stage = 0; stage < this->stepper.get_num_stages(); ++stage) {
            this->mesh.CallForEachElement(volume_kernel);

            this->mesh.CallForEachElement(source_kernel);

            this->mesh.CallForEachInterface(interface_kernel);

            this->mesh.CallForEachBoundary(boundary_kernel);

            this->mesh.CallForEachElement(update_kernel);

            this->mesh.CallForEachElement(scrutinize_solution_kernel);

            ++(this->stepper);
        }

        this->mesh.CallForEachElement(swap_states_kernel);

        if (step % 360 == 0) {
            log_file << "Step: " << this->stepper.get_step() << std::endl;

            ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
            ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
        }
    }

    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    mesh.CallForEachElement(compute_residual_L2_kernel);

    log_file << "residual L2 norm: " << sqrt(residual_L2) << std::endl;
}

#endif
