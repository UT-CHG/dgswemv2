#include "swe_parser.hpp"

namespace SWE {
Parser::Parser(const InputParameters<SWE::Inputs>& input) {
    if (input.problem_input.meteo_forcing.type == MeteoForcingType::Enable) {
        this->parsing_input = true;

        this->meteo_parse_frequency =
            (uint)std::ceil(input.problem_input.meteo_forcing.frequency / input.stepper_input.dt);
        this->meteo_data_file = input.problem_input.meteo_forcing.meteo_data_file;
    }
}

Parser::Parser(const InputParameters<SWE::Inputs>& input, const uint locality_id, const uint submesh_id)
    : Parser(input) {}  // this is for partitioned input files

void Parser::ParseMeteoInput(const Stepper& stepper) {
    uint step = stepper.GetStep();

    if (this->node_meteo_data_step.find(step - this->meteo_parse_frequency) != this->node_meteo_data_step.end()) {
        this->node_meteo_data_step.erase(this->node_meteo_data_step.find(step - this->meteo_parse_frequency));
    }

    if (this->node_meteo_data_step.find(step) == this->node_meteo_data_step.end()) {
        std::string meteo_data_file_name = this->meteo_data_file;

        meteo_data_file_name.insert(meteo_data_file_name.find_last_of("."), '_' + std::to_string(step));

        if (!Utilities::file_exists(meteo_data_file_name)) {
            throw std::logic_error("Fatal Error: meteo data file " + meteo_data_file_name + " was not found!\n");
        }

        std::ifstream meteo_file(meteo_data_file_name);

        uint node_id;
        std::vector<double> meteo_data(3);

        std::string line;
        while (std::getline(meteo_file, line)) {
            std::istringstream input_string(line);

            if (!(input_string >> node_id >> meteo_data[0] >> meteo_data[1] >> meteo_data[2]))
                break;

            this->node_meteo_data_step[step][node_id] = meteo_data;
        }
    }

    if (this->node_meteo_data_step.find(step + this->meteo_parse_frequency) == this->node_meteo_data_step.end()) {
        std::string meteo_data_file_name = this->meteo_data_file;

        meteo_data_file_name.insert(meteo_data_file_name.find_last_of("."),
                                    '_' + std::to_string(step + this->meteo_parse_frequency));

        if (!Utilities::file_exists(meteo_data_file_name)) {
            throw std::logic_error("Fatal Error: meteo data file " + meteo_data_file_name + " was not found!\n");
        }

        std::ifstream meteo_file(meteo_data_file_name);

        uint node_id;
        std::vector<double> meteo_data(3);

        std::string line;
        while (std::getline(meteo_file, line)) {
            std::istringstream input_string(line);

            if (!(input_string >> node_id >> meteo_data[0] >> meteo_data[1] >> meteo_data[2]))
                break;

            this->node_meteo_data_step[step + this->meteo_parse_frequency][node_id] = meteo_data;
        }
    }
}

void Parser::InterpolateMeteoData(const Stepper& stepper) {
    uint step = stepper.GetStep();

    uint step_start = step - step % this->meteo_parse_frequency;
    uint step_end   = step_start + this->meteo_parse_frequency;

    double t_start = step_start * stepper.GetDT();
    double t_end   = step_end * stepper.GetDT();

    double interp_factor = (stepper.GetTimeAtCurrentStage() - t_start) / (t_end - t_start);

    for (auto it = this->node_meteo_data_interp.begin(); it != this->node_meteo_data_interp.end(); ++it) {
        it->second[0] = this->node_meteo_data_step[step_start][it->first][0] +
                        interp_factor * (this->node_meteo_data_step[step_end][it->first][0] -
                                         this->node_meteo_data_step[step_start][it->first][0]);

        it->second[1] = this->node_meteo_data_step[step_start][it->first][1] +
                        interp_factor * (this->node_meteo_data_step[step_end][it->first][1] -
                                         this->node_meteo_data_step[step_start][it->first][1]);

        it->second[2] = this->node_meteo_data_step[step_start][it->first][2] +
                        interp_factor * (this->node_meteo_data_step[step_end][it->first][2] -
                                         this->node_meteo_data_step[step_start][it->first][2]);
    }
}
}