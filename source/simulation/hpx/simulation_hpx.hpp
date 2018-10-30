#ifndef SIMULATION_HPX_HPP
#define SIMULATION_HPX_HPP

#include "general_definitions.hpp"
#include "preprocessor/input_parameters.hpp"
#include "utilities/file_exists.hpp"
#include "sim_unit_hpx_base.hpp"

#include <hpx/util/unwrapped.hpp>
//#include "simulation/hpx/load_balancer/load_balancer_headers.hpp"

template <typename ClientType>
hpx::future<double> ComputeL2Residual(std::vector<ClientType>& clients) {
    std::vector<hpx::future<double>> res_futures;

    for (uint id = 0; id < clients.size(); id++) {
        res_futures.push_back(clients[id].ResidualL2());
    }

    return hpx::when_all(res_futures).then([](auto&& res_futures) -> double {
        std::vector<double> res = hpx::util::unwrap(res_futures.get());
        double combined_res{0};
        for (double r : res) {
            combined_res += r;
        }
        return combined_res;
    });
}

class HPXSimulation : public hpx::components::component_base<HPXSimulation> {
  public:
    using ClientType = HPXSimulationUnitClient;

  private:
    uint n_steps;

    std::vector<ClientType> simulation_unit_clients;

  public:
    HPXSimulation() = default;
    HPXSimulation(const std::string& input_string);

    hpx::future<void> Run();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, Run, RunAction);

    hpx::future<double> ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulation, ResidualL2, ResidualL2Action);
};

HPXSimulation::HPXSimulation(const std::string& input_string) {
    const uint locality_id          = hpx::get_locality_id();
    const hpx::naming::id_type here = hpx::find_here();

    InputParameters<> input(input_string);

    this->n_steps = (uint)std::ceil(input.stepper_input.run_time / input.stepper_input.dt);

    hpx::future<void> lb_future = hpx::make_ready_future();
//        LoadBalancer::AbstractFactory::initialize_locality_and_world_models<ProblemType>(locality_id, input_string);

    std::string submesh_file_prefix =
        input.mesh_input.mesh_file_name.substr(0, input.mesh_input.mesh_file_name.find_last_of('.')) + "_" +
        std::to_string(locality_id) + '_';
    std::string submesh_file_postfix = input.mesh_input.mesh_file_name.substr(
        input.mesh_input.mesh_file_name.find_last_of('.'), input.mesh_input.mesh_file_name.size());

    std::vector<hpx::future<void>> registration_futures;

    uint submesh_id = 0;
    while (Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix)) {
        this->simulation_unit_clients.emplace_back( HPXSimulationUnitFactory::Create( here, input_string, locality_id, submesh_id)
            );

        registration_futures.push_back(this->simulation_unit_clients.back().register_as(
            std::string{ClientType::GetBasename()} + std::to_string(locality_id) + '_' + std::to_string(submesh_id)));

        ++submesh_id;
    }

    hpx::when_all(registration_futures).get();

    lb_future.get();
}

hpx::future<void> HPXSimulation::Run() {
    std::vector<hpx::future<void>> simulation_futures;

    for (auto& sim_unit_client : this->simulation_unit_clients) {
        simulation_futures.push_back(sim_unit_client.Preprocessor());
    }

    for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
        simulation_futures[sim_id] = simulation_futures[sim_id].then(
            [this, sim_id](auto&&) { return this->simulation_unit_clients[sim_id].Launch(); });
    }

    for (uint step = 1; step <= this->n_steps; step++) {
        for (uint sim_id = 0; sim_id < this->simulation_unit_clients.size(); sim_id++) {
            simulation_futures[sim_id] = simulation_futures[sim_id].then([this, sim_id](auto&& f) {
                f.get();  // check for exceptions
                return this->simulation_unit_clients[sim_id].Step();
            });
        }
    }

    return hpx::when_all(simulation_futures); /*.then([](auto&&) {
        LoadBalancer::AbstractFactory::reset_locality_and_world_models<ProblemType>();
        });*/
}

hpx::future<double> HPXSimulation::ResidualL2() {
    return ComputeL2Residual(this->simulation_unit_clients);
}

class HPXSimulationClient : hpx::components::client_base<HPXSimulationClient, HPXSimulation> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationClient, HPXSimulation>;

  public:
    HPXSimulationClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}
    HPXSimulationClient(hpx::id_type&& id) : BaseType(std::move(id)) {}

    hpx::future<void> Run() {
        using ActionType = typename HPXSimulation::RunAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<double> ResidualL2() {
        using ActionType = typename HPXSimulation::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }
};

using hpx_simulation_swe_component_ = hpx::components::simple_component<HPXSimulation>;
HPX_REGISTER_COMPONENT(hpx_simulation_swe_component_, hpx_simulation_swe_);

#endif
