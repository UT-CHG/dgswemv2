#ifndef SIMULATION_OMPI_HPP
#define SIMULATION_OMPI_HPP

#include <omp.h>

#include "general_definitions.hpp"
#include "preprocessor/input_parameters.hpp"
#include "utilities/file_exists.hpp"
#include "sim_unit_ompi.hpp"

template <typename ProblemType>
class OMPISimulation : public OMPISimulationBase {
  private:
    uint n_steps;

    std::vector<std::unique_ptr<OMPISimulationUnit<ProblemType>>> sim_units;
    typename ProblemType::ProblemGlobalDataType global_data;

    typename ProblemType::ProblemStepperType stepper;

  public:
    OMPISimulation() = default;
    OMPISimulation(const std::string& input_string);

    void Run();
    void ComputeL2Residual();
    void Finalize();
};

template <typename ProblemType>
OMPISimulation<ProblemType>::OMPISimulation(const std::string& input_string) {
    int locality_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

    InputParameters<typename ProblemType::ProblemInputType> input(input_string);

    this->n_steps = (uint)std::ceil(input.stepper_input.run_time / input.stepper_input.dt);

    this->stepper = typename ProblemType::ProblemStepperType(input.stepper_input);

    std::string submesh_file_prefix =
        input.mesh_input.mesh_file_name.substr(0, input.mesh_input.mesh_file_name.find_last_of('.')) + "_" +
        std::to_string(locality_id) + '_';
    std::string submesh_file_postfix = input.mesh_input.mesh_file_name.substr(
        input.mesh_input.mesh_file_name.find_last_of('.'), input.mesh_input.mesh_file_name.size());

    uint submesh_id = 0;

    while (Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix)) {
        this->sim_units.emplace_back(new OMPISimulationUnit<ProblemType>(input_string, locality_id, submesh_id));

        ++submesh_id;
    }

    if (this->sim_units.empty()) {
        std::cerr << "Warning: MPI Rank " << locality_id << " has not been assigned any work. This may inidicate\n"
                  << "         poor partitioning and imply degraded performance." << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulation<ProblemType>::Run() {
#pragma omp parallel
    {
        uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

        n_threads = (uint)omp_get_num_threads();
        thread_id = (uint)omp_get_thread_num();

        sim_per_thread = (this->sim_units.size() + n_threads - 1) / n_threads;

        begin_sim_id = sim_per_thread * thread_id;
        end_sim_id   = std::min(sim_per_thread * (thread_id + 1), (uint)this->sim_units.size());

        ProblemType::preprocessor_ompi(this->sim_units, this->global_data, this->stepper, begin_sim_id, end_sim_id);

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (this->sim_units[su_id]->writer.WritingLog()) {
                this->sim_units[su_id]->writer.GetLogFile() << std::endl
                                                            << "Launching Simulation!" << std::endl
                                                            << std::endl;
            }

            if (this->sim_units[su_id]->writer.WritingOutput()) {
                this->sim_units[su_id]->writer.WriteFirstStep(this->stepper,
                                                              this->sim_units[su_id]->discretization.mesh);
            }
        }

        for (uint step = 1; step <= this->n_steps; ++step) {
            ProblemType::step_ompi(this->sim_units, this->global_data, this->stepper, begin_sim_id, end_sim_id);
        }
    }  // close omp parallel region
}

template <typename ProblemType>
void OMPISimulation<ProblemType>::ComputeL2Residual() {
    int locality_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

    double global_l2{0};
    double residual_l2{0};

    for (auto& sim_unit : this->sim_units) {
        sim_unit->discretization.mesh.CallForEachElement([this, &sim_unit, &residual_l2](auto& elt) {
            residual_l2 += ProblemType::compute_residual_L2(this->stepper, elt);
        });
    }

    MPI_Reduce(&residual_l2, &global_l2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (locality_id == 0) {
        std::cout << "L2 error: " << std::setprecision(15) << std::sqrt(global_l2) << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulation<ProblemType>::Finalize() {
    ProblemType::finalize_simulation(this->global_data);
}

#endif