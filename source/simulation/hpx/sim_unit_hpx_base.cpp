#include "sim_unit_hpx_base.hpp"
#include "sim_unit_hpx.hpp"

#include <yaml-cpp/yaml.h>

HPXSimulationUnitClient HPXSimulationUnitFactory::Create(const hpx::naming::id_type& here,
                                                         const std::string& input_string,
                                                         const uint locality_id,
                                                         const uint submesh_id) {
    YAML::Node input_ = YAML::LoadFile(input_string);

    if ( input_["problem"] && input_["problem"]["name"]) {
        std::string problem_name = input_["problem"]["name"].as<std::string>();

        if ( problem_name == "rkdg_swe" ) {
            return hpx::components::new_<HPXSimulationUnit<SWE::RKDG::Problem>>(here, input_string, locality_id, submesh_id);
        } else {
            throw std::runtime_error{"Unknown problem name: "+problem_name};
        }
    } else {
        throw std::runtime_error{"Input file has no node 'problem' with field 'name'\n"};
    }
}