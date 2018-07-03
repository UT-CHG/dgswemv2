#ifndef IHDG_SIMULATION_OMPI_HPP
#define IHDG_SIMULATION_OMPI_HPP

#include "general_definitions.hpp"

#include "preprocessor/input_parameters.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"
#include "communication/ompi_communicator.hpp"
#include "utilities/file_exists.hpp"

#include "ihdg_sim_unit_ompi.hpp"
#include "simulation/writer.hpp"

namespace IHDG {
template <typename ProblemType>
class OMPISimulation {
  private:
    uint n_steps;
    uint n_stages;

    std::vector<std::unique_ptr<OMPISimulationUnit<ProblemType>>> simulation_units;

  public:
    OMPISimulation() = default;
    OMPISimulation(const std::string& input_string);

    void Run();
    void ComputeL2Residual();
};

template <typename ProblemType>
OMPISimulation<ProblemType>::OMPISimulation(const std::string& input_string) {
    int locality_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

    InputParameters<typename ProblemType::ProblemInputType> input(input_string);

    this->n_steps  = (uint)std::ceil(input.stepper_input.run_time / input.stepper_input.dt);
    this->n_stages = input.stepper_input.nstages;

    std::string submesh_file_prefix =
        input.mesh_input.mesh_file_name.substr(0, input.mesh_input.mesh_file_name.find_last_of('.')) + "_" +
        std::to_string(locality_id) + '_';
    std::string submesh_file_postfix = input.mesh_input.mesh_file_name.substr(
        input.mesh_input.mesh_file_name.find_last_of('.'), input.mesh_input.mesh_file_name.size());

    uint submesh_id = 0;

    while (Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix)) {
        this->simulation_units.emplace_back(new OMPISimulationUnit<ProblemType>(input_string, locality_id, submesh_id));

        ++submesh_id;
    }
}

template <typename ProblemType>
void OMPISimulation<ProblemType>::Run() {
#pragma omp parallel
    {
        uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

        n_threads = (uint)omp_get_num_threads();
        thread_id = (uint)omp_get_thread_num();

        sim_per_thread = (this->simulation_units.size() + n_threads - 1) / n_threads;

        begin_sim_id = sim_per_thread * thread_id;
        end_sim_id   = std::min(sim_per_thread * (thread_id + 1), (uint)this->simulation_units.size());

        for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
            this->simulation_units[sim_unit_id]->PostReceivePrerocStage();
        }

        for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
            this->simulation_units[sim_unit_id]->WaitAllPreprocSends();
        }

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

                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->WaitAllSends();
                }
            }

            for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                this->simulation_units[sim_unit_id]->Step();
            }
        }
    }  // close omp parallel region
}

template <typename ProblemType>
void OMPISimulation<ProblemType>::ComputeL2Residual() {
    int locality_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

    double global_l2{0};

    double residual_l2{0};
    for (auto& sim_unit : this->simulation_units) {
        residual_l2 += sim_unit->ResidualL2();
    }

    MPI_Reduce(&residual_l2, &global_l2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (locality_id == 0) {
        std::cout << "L2 error: " << std::setprecision(14) << std::sqrt(global_l2) << std::endl;
    }
}
}

#endif
