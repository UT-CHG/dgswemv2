#ifndef GN_PARSER_HPP
#define GN_PARSER_HPP

#include "preprocessor/input_parameters.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"

#include "problem/Green-Naghdi/gn_definitions.hpp"

#include "utilities/file_exists.hpp"

namespace GN {
class Parser : public SWE::Parser {
  public:
    Parser() = default;
    Parser(const InputParameters<GN::Inputs>& input);
    Parser(const InputParameters<GN::Inputs>& input, const uint locality_id, const uint submesh_id);

    template <typename MeshType>
    void ParseInput(const RKStepper& stepper, MeshType& mesh);

  public:
#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & parsing_input
            & meteo_parse_frequency
            & meteo_data_file
            & node_meteo_data_step
            & node_meteo_data_interp;
        // clang-format on
    }
#endif
};

template <typename MeshType>
void Parser::ParseInput(const RKStepper& stepper, MeshType& mesh) {
    if (SWE::SourceTerms::meteo_forcing) {
        if (stepper.GetStep() % this->meteo_parse_frequency == 0 && stepper.GetStage() == 0) {
            this->ParseMeteoInput(stepper);
        }

        // Initialize container to store parsed data and store pointers for fast access
        if (stepper.GetStep() == 0 && stepper.GetStage() == 0) {
            mesh.CallForEachElement([this](auto& elt) {
                const std::vector<uint>& node_ID = elt.GetNodeID();

                for (uint node = 0; node < elt.data.get_nnode(); ++node) {
                    if (this->node_meteo_data_interp.find(node_ID[node]) == this->node_meteo_data_interp.end()) {
                        this->node_meteo_data_interp[node_ID[node]] = std::vector<double>(3);
                    }

                    elt.data.source.parsed_meteo_data[node] = &this->node_meteo_data_interp[node_ID[node]];
                }
            });
        }

        this->InterpolateMeteoData(stepper);

        mesh.CallForEachElement([this](auto& elt) {
            for (uint node = 0; node < elt.data.get_nnode(); ++node) {
                elt.data.source.tau_s[node][GlobalCoord::x] = elt.data.source.parsed_meteo_data[node]->at(0);
                elt.data.source.tau_s[node][GlobalCoord::y] = elt.data.source.parsed_meteo_data[node]->at(1);
                elt.data.source.p_atm[node]                 = elt.data.source.parsed_meteo_data[node]->at(2);
            }
        });
    }
}
}

#endif