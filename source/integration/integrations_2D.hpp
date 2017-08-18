#ifndef INTEGRATIONS_2D_HPP
#define INTEGRATIONS_2D_HPP

#include "../general_definitions.hpp"

namespace Integration {
class Dunavant_2D : Integration<2> {
  public:
    std::pair<std::vector<double>, std::vector<Point<2>>> GetRule(uint);
    uint GetNumGP(uint p);

  private:
    std::vector<uint> PermutationData(uint);
    std::pair<std::vector<double>, std::vector<Point<3>>> GPData(uint);
};
}

#endif