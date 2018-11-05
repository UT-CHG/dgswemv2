#include "sim_unit_hpx_base.hpp"
#include "sim_unit_hpx.hpp"

#include <yaml-cpp/yaml.h>

HPX_REGISTER_COMPONENT_HEAP(hpx::components::managed_component<HPXSimulationUnitBase>);
HPX_DEFINE_GET_COMPONENT_TYPE(HPXSimulationUnitBase);

HPXSimulationUnitClient HPXSimulationUnitFactory::Create(const hpx::naming::id_type& here,
                                                         const std::string& input_string,
                                                         const uint locality_id,
                                                         const uint submesh_id) {
    YAML::Node input_ = YAML::LoadFile(input_string);

    if ( input_["problem"] && input_["problem"]["name"]) {
        std::string problem_name = input_["problem"]["name"].as<std::string>();

        if ( problem_name == "rkdg_swe" ) {
            return HPXSimulationUnitFactory::CreateSimulationUnit<SWE::RKDG::Problem>(here,input_string, locality_id, submesh_id);
        } else if ( problem_name == "ehdg_swe" ) {
            return HPXSimulationUnitFactory::CreateSimulationUnit<SWE::EHDG::Problem>(here,input_string, locality_id, submesh_id);
        } else {
            throw std::runtime_error{"Unknown problem name: "+problem_name};
        }
    } else {
        throw std::runtime_error{"Input file has no node 'problem' with field 'name'\n"};
    }
}