#ifndef SWE_PARSER_HPP
#define SWE_PARSER_HPP

#include "utilities/file_exists.hpp"
#include "preprocessor/input_parameters.hpp"

namespace SWE {
class Parser {
  private:
    bool parsing_input = false;

    uint meteo_parse_frequency;
    std::string meteo_data_file;
    std::map<uint, std::map<uint, std::vector<double>>> node_meteo_data_step;
    std::map<uint, std::vector<double>> node_meteo_data_interp;

  public:
    Parser() = default;
    template <typename ProblemSpecificInputType>
    Parser(const InputParameters<ProblemSpecificInputType>& input);
    template <typename ProblemSpecificInputType>
    Parser(const InputParameters<ProblemSpecificInputType>& input, const uint locality_id, const uint submesh_id);

    bool ParsingInput() { return parsing_input; }

    template <typename StepperType, typename MeshType>
    void ParseInput(const StepperType& stepper, MeshType& mesh);

  private:
    template <typename StepperType>
    void ParseMeteoInput(const StepperType& stepper);
    template <typename StepperType>
    void InterpolateMeteoData(const StepperType& stepper);

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

template <typename ProblemSpecificInputType>
Parser::Parser(const InputParameters<ProblemSpecificInputType>& input) {
    if (input.problem_input.meteo_forcing.type == MeteoForcingType::Enable) {
        this->parsing_input = true;

        this->meteo_parse_frequency =
            (uint)std::ceil(input.problem_input.meteo_forcing.frequency / input.stepper_input.dt);
        this->meteo_data_file = input.problem_input.meteo_forcing.meteo_data_file;
    }
}

template <typename ProblemSpecificInputType>
Parser::Parser(const InputParameters<ProblemSpecificInputType>& input, const uint locality_id, const uint submesh_id)
    : Parser(input) {}  // this is for partitioned input files

template <typename StepperType, typename MeshType>
void Parser::ParseInput(const StepperType& stepper, MeshType& mesh) {
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

template <typename StepperType>
void Parser::ParseMeteoInput(const StepperType& stepper) {
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

template <typename StepperType>
void Parser::InterpolateMeteoData(const StepperType& stepper) {
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

#endif