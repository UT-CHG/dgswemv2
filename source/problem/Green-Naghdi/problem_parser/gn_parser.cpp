#include "gn_parser.hpp"

namespace GN {
Parser::Parser(const InputParameters<GN::Inputs>& input) {
    if (input.problem_input.meteo_forcing.type == SWE::MeteoForcingType::Enable) {
        this->parsing_input = true;

        this->meteo_parse_frequency =
            (uint)std::ceil(input.problem_input.meteo_forcing.frequency / input.stepper_input.dt);
        this->meteo_data_file = input.problem_input.meteo_forcing.meteo_data_file;
    }
}

Parser::Parser(const InputParameters<GN::Inputs>& input, const uint locality_id, const uint submesh_id)
    : Parser(input) {}  // this is for partitioned input files
}