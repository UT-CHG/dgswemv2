#ifndef SERIAL_SIMULATION_HPP
#define SERIAL_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"

#include "simulation.hpp"

template <typename ProblemType>
class SerialSimulation : public Simulation<ProblemType> {
public:
  SerialSimulation() : Simulation<ProblemType>() {}
  SerialSimulation(std::string& input_string) : Simulation<ProblemType>(input_string) {}

  void Run();
};

template <typename ProblemType>
void SerialSimulation<ProblemType>::Run() {
    // we write these gross looking wrapper functions to append the stepper in a
    // way that allows us to keep the
    // the nice std::for_each notation without having to define stepper within
    // each element
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
            this->log_file << "Step: " << this->stepper.get_step() << std::endl;
            ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
            ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
        }
    }
}

#endif
