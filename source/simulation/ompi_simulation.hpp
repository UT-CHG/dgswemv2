#ifndef OMPI_SIMULATION_HPP
#define OMPI_SIMULATION_HPP

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "../communication/ompi_communicator.hpp"

#include "utilities/file_exists.hpp"

template <typename ProblemType>
class OMPISimulationUnit {
  private:
    InputParameters input;

    Stepper stepper;
    typename ProblemType::ProblemMeshType mesh;
    OMPICommunicator communicator;

    std::string log_file_name;

  public:
    OMPISimulationUnit() : input(), stepper(input.rk.nstages, input.rk.order, input.dt), mesh(input.polynomial_order) {}
    OMPISimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id)
        : input(input_string, locality_id, submesh_id),
          stepper(input.rk.nstages, input.rk.order, input.dt),
          mesh(input.polynomial_order),
          communicator(std::string(input.mesh_file_name + "_meta"), locality_id, submesh_id) {
        input.ReadMesh();

        mesh.SetMeshName(input.mesh_data.mesh_name);

        this->log_file_name = "output/" + input.mesh_data.mesh_name + "_log";

        std::ofstream log_file(this->log_file_name, std::ofstream::out);

        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << input.mesh_file_name
                 << " mesh" << std::endl << std::endl;

        initialize_mesh<ProblemType, OMPICommunicator>(this->mesh, input.mesh_data, communicator);
    }

    void Launch();
    void ExchangeData();
    void PreReceiveStage();
    void PostReceiveStage();
    void Step();
};

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::Launch() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);

    log_file << std::endl << "Launching Simulation!" << std::endl << std::endl;

    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);

    ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
    ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::ExchangeData() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
    log_file << "Current (time, stage): (" << this->stepper.get_t_at_curr_stage() << ',' << this->stepper.get_stage()
             << ')' << std::endl;
#endif
    auto distributed_boundary_send_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_send_kernel(this->stepper, dbound);
    };

    this->communicator.ReceiveAll(this->stepper.get_timestamp());

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_send_kernel);

    this->communicator.SendAll(this->stepper.get_timestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PreReceiveStage() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
    log_file << "Starting work before receive" << std::endl;
#endif
    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    this->mesh.CallForEachElement(volume_kernel);

    this->mesh.CallForEachElement(source_kernel);

    this->mesh.CallForEachInterface(interface_kernel);
#ifdef VERBOSE
    log_file << "Finished work before receive" << std::endl;
#endif
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PostReceiveStage() {
    std::ofstream log_file(this->log_file_name, std::ofstream::app);
#ifdef VERBOSE
    log_file << "Starting to wait on receive with timestamp: " << this->stepper.get_timestamp() << std::endl;
#endif
    this->communicator.WaitAllReceives(this->stepper.get_timestamp());
#ifdef VERBOSE
    log_file << "Starting work after receive" << std::endl;
#endif
    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    auto distributed_boundary_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_kernel(this->stepper, dbound);
    };

    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

    auto scrutinize_solution_kernel = [this](auto& elt) {
        ProblemType::scrutinize_solution_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachBoundary(boundary_kernel);

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_kernel);

    this->mesh.CallForEachElement(update_kernel);

    this->mesh.CallForEachElement(scrutinize_solution_kernel);

    ++(this->stepper);
#ifdef VERBOSE
    log_file << "Finished work after receive" << std::endl << std::endl;
#endif
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::Step() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->stepper.get_step() % 360 == 0) {
        std::ofstream log_file(this->log_file_name, std::ofstream::app);

        log_file << "Step: " << this->stepper.get_step() << std::endl;

        ProblemType::write_VTK_data_kernel(this->stepper, this->mesh);
        ProblemType::write_modal_data_kernel(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
class OMPISimulation {
  private:
    uint n_steps;
    uint n_stages;

    std::vector<OMPISimulationUnit<ProblemType>*> simulation_units;

  public:
    OMPISimulation() = default;
    OMPISimulation(const std::string& input_string) {
        int locality_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

        InputParameters input(input_string);

        this->n_steps = (uint)std::ceil(input.T_end / input.dt);
        this->n_stages = input.rk.nstages;

        std::string mesh_file_prefix = input.mesh_file_name;
        mesh_file_prefix.erase(mesh_file_prefix.size() - 3);
        mesh_file_prefix += '_' + std::to_string(locality_id) + '_';

        uint submesh_id = 0;
        while (Utilities::file_exists(mesh_file_prefix + std::to_string(submesh_id) + ".14")) {
            this->simulation_units.push_back(
                new OMPISimulationUnit<ProblemType>(input_string, locality_id, submesh_id));

            ++submesh_id;
        }
    }

    void Run();
};

template <typename ProblemType>
void OMPISimulation<ProblemType>::Run() {
#pragma omp parallel
    {
        uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

        n_threads = (uint)omp_get_num_threads();
        thread_id = (uint)omp_get_thread_num();

        sim_per_thread = this->simulation_units.size() / n_threads;

        begin_sim_id = sim_per_thread * thread_id;
        end_sim_id = sim_per_thread * (thread_id + 1);

        if (end_sim_id > this->simulation_units.size())
            end_sim_id = this->simulation_units.size();

        for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
            this->simulation_units[sim_unit_id]->Launch();
        }

        for (uint step = 1; step <= this->n_steps; step++) {
            for (uint stage = 0; stage < this->n_stages; stage++) {
                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->ExchangeData();
                }

                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->PreReceiveStage();
                }

                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->PostReceiveStage();
                }
            }

            for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                this->simulation_units[sim_unit_id]->Step();
            }
        }
    }
}

#endif
