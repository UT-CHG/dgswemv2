#include "simulation/serial/simulation_base.hpp"
#include "simulation/serial/simulation.hpp"

namespace Serial {

std::unique_ptr<SimulationBase> SimulationFactory::Create(const std::string& input_string) {
    YAML::Node input_ = YAML::LoadFile(input_string);

    if (input_["problem"] && input_["problem"]["name"] ) {
        std::string problem_name = input_["problem"]["name"].as<std::string>();

        if ( problem_name == "rkdg_swe") {
            return SimulationFactory::CreateSimulation<SWE::RKDG::Problem>(input_string);
        } else if (problem_name == "ehdg_swe") {
            return SimulationFactory::CreateSimulation<SWE::EHDG::Problem>(input_string);
        } else if (problem_name == "ihdg_swe") {
            return SimulationFactory::CreateSimulation<SWE::IHDG::Problem>(input_string);
        } else {
            throw std::runtime_error{"Unknown problem name: "+problem_name};
        }
    } else {
        throw std::runtime_error{"Input file has no node 'problem with field 'name'\n"};
    }
}

}