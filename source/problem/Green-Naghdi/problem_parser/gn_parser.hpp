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
    Parser(const InputParameters<GN::Inputs>& input) : SWE::Parser(input) {}
    Parser(const InputParameters<GN::Inputs>& input, const uint locality_id, const uint submesh_id)
        : SWE::Parser(input, locality_id, submesh_id) {}
};
}

#endif