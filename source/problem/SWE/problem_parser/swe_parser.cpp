#include "swe_parser.hpp"

namespace SWE {
void Parser::ParseMeteoInput(uint step) {
    if (this->meteo_forcing_type == SWE::MeteoForcingType::Test) {
        if (this->node_meteo_data_step.find(step - this->meteo_parse_frequency) != this->node_meteo_data_step.end()) {
            this->node_meteo_data_step.erase(this->node_meteo_data_step.find(step - this->meteo_parse_frequency));
        }

        if (this->node_meteo_data_step.find(step) == this->node_meteo_data_step.end()) {
            std::string meteo_data_file_name = this->meteo_data_file;

            meteo_data_file_name.insert(meteo_data_file_name.find_last_of("."), '_' + std::to_string(step));

            std::ifstream meteo_file(meteo_data_file_name);

            if (!meteo_file) {
                std::string err_msg = "Fatal Error: Meteo data file " + meteo_data_file_name + " was not found\n";
                throw std::logic_error(err_msg);
            }

            uint                node_id;
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

            std::ifstream meteo_file(meteo_data_file_name);

            if (!meteo_file) {
                std::string err_msg = "Fatal Error: Meteo data file " + meteo_data_file_name + " not found\n";
                throw std::logic_error(err_msg);
            }

            uint                node_id;
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
}

void Parser::CalculateMeteoData(uint step) {
    if (this->meteo_forcing_type == SWE::MeteoForcingType::Test) {
        uint step_begin = step - step % this->meteo_parse_frequency;
        uint step_end   = step_begin + this->meteo_parse_frequency;

        double interp_factor = step % this->meteo_parse_frequency / ((double)this->meteo_parse_frequency);

        this->node_meteo_data = this->node_meteo_data_step[step_begin];

        for (auto it = this->node_meteo_data.begin(); it != this->node_meteo_data.end(); ++it) {
            it->second[0] += interp_factor * (this->node_meteo_data_step[step_end][it->first][0] -
                                              this->node_meteo_data_step[step_begin][it->first][0]);

            it->second[1] += interp_factor * (this->node_meteo_data_step[step_end][it->first][1] -
                                              this->node_meteo_data_step[step_begin][it->first][1]);

            it->second[2] += interp_factor * (this->node_meteo_data_step[step_end][it->first][2] -
                                              this->node_meteo_data_step[step_begin][it->first][2]);
        }
    }
}
}