#ifndef ADCIRC_FORMAT_HPP
#define ADCIRC_FORMAT_HPP

#include "../../general_definitions.hpp"
#include "../../problem/SWE/swe_definitions.hpp"

#include "../../shape/shapes_2D.hpp"

struct AdcircFormat {
    std::string name;
    std::unordered_map<int, std::array<double, 3>> nodes;
    std::unordered_map<int, std::array<int, 4>> elements;

    // see
    // http://adcirc.org/home/documentation/users-manual-v51/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14/
    int NOPE;  // number of open boundaries
    int NETA;  // total number of nodes for open boundary

    std::vector<std::vector<int>> NBDV;  // Ordered node numbers

    int NBOU;
    int NVEL;

    std::vector<int> IBTYPE;             // boundary type
    std::vector<std::vector<int>> NBVV;  // node numbers on normal flow boundary
                                         // segment k

    AdcircFormat(const std::string& in_name);

    void write_to(const char* out_name) const;

    SWE::BoundaryConditions get_ibtype(std::array<int, 2>& node_pair) const;

    bool has_edge(std::vector<int>::const_iterator cbegin,
                  std::vector<int>::const_iterator cend,
                  std::array<int, 2>& node_pair) const;
};

#endif
