#ifndef SWE_PARSER_HPP
#define SWE_PARSER_HPP

#include "preprocessor/input_parameters.hpp"

namespace SWE {
class Parser {
  private:
    bool parsing_input;

    SWE::MeteoForcingType meteo_forcing_type;
    uint meteo_parse_frequency;
    std::string meteo_data_file;
    std::map<uint, std::vector<double>> node_meteo_data;

  public:
    Parser() = default;
    Parser(const InputParameters<SWE::Inputs>& input);
    Parser(const InputParameters<SWE::Inputs>& input,
           const uint locality_id,
           const uint submesh_id);  // this is for partitioned input files

    bool ParsingInput() { return parsing_input; }

    template <typename MeshType>
    void ParseInput(const Stepper& stepper, MeshType& mesh);

  private:
    void ParseMeteoInput();

#ifdef HAS_HPX
  public:
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

Parser::Parser(const InputParameters<SWE::Inputs>& input) {
    this->parsing_input = input.problem_input.parse_input;

    meteo_forcing_type = input.problem_input.meteo_forcing.type;
    meteo_parse_frequency = (uint)std::ceil(input.problem_input.meteo_forcing.frequency / input.dt);
    this->meteo_data_file = input.problem_input.meteo_forcing.meteo_data_file;
}

Parser::Parser(const InputParameters<SWE::Inputs>& input, const uint locality_id, const uint submesh_id)
    : Parser(input) {}

template <typename MeshType>
void Parser::ParseInput(const Stepper& stepper, MeshType& mesh) {
    if (SWE::SourceTerms::meteo_forcing) {
        if (stepper.get_step() % this->meteo_parse_frequency == 0) {
            this->ParseMeteoInput();
        }

        mesh.CallForEachElement([this](auto& elt) {
            std::vector<uint>& node_ID = elt.GetNodeID();

            //# of node != # of vrtx in case we have an iso-p element with p>1
            // I assume we will have values only at vrtx in files
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                elt.data.source.tau_s[GlobalCoord::x][vrtx] = this->node_meteo_data[node_ID[vrtx]][0];
                elt.data.source.tau_s[GlobalCoord::y][vrtx] = this->node_meteo_data[node_ID[vrtx]][1];

                elt.data.source.p_atm[vrtx] = this->node_meteo_data[node_ID[vrtx]][2];
            }
        });
    }
}

void Parser::ParseMeteoInput() {
    if (this->meteo_forcing_type == SWE::MeteoForcingType::Test) {
        std::ifstream meteo_file(this->meteo_data_file);

        uint node_id;
        std::vector<double> meteo_data(3);

        std::string line;
        while (std::getline(meteo_file, line)) {
            std::istringstream input_string(line);

            if (!(input_string >> node_id >> meteo_data[0] >> meteo_data[1] >> meteo_data[2]))
                break;

            this->node_meteo_data[node_id] = meteo_data;
        }
    }
}

#ifdef HAS_HPX
template <typename Archive>
void Parser::serialize(Archive& ar, unsigned) {
    ar & meteo_forcing_type
       & meteo_parse_frequency
       & meteo_data_file
       & node_meteo_data
       & parsing_input;
}
#endif
}
#endif