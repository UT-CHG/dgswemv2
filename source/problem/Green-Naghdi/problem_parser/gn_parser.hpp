#ifndef GN_PARSER_HPP
#define GN_PARSER_HPP

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