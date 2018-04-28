#ifndef ADCIRC_FORMAT_HPP
#define ADCIRC_FORMAT_HPP

#include "../../general_definitions.hpp"
#include "../../problem/SWE/swe_definitions.hpp"

#include "../../shape/shapes_2D.hpp"

struct AdcircFormat {
    std::string                                     name;
    std::unordered_map<uint, std::array<double, 3>> nodes;
    std::unordered_map<uint, std::array<uint, 4>>   elements;

    // see
    // http://adcirc.org/home/documentation/users-manual-v51/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14/
    uint NOPE;  // number of open boundaries
    uint NETA;  // total number of nodes for open boundary

    std::vector<std::vector<uint>> NBDV;  // Ordered node numbers

    uint NBOU;
    uint NVEL;

    std::vector<uint>              IBTYPE;  // boundary type
    std::vector<std::vector<uint>> NBVV;    // node numbers on normal flow boundary segment k

    std::map<uint, std::vector<uint>>   IBCONN;     // node numbers of back nodes for internal barrier segment k
    std::map<uint, std::vector<double>> BARINTH;    // for internal barrier segment k
    std::map<uint, std::vector<double>> BARINCFSB;  // for internal barrier segment k
    std::map<uint, std::vector<double>> BARINCFSP;  // for internal barrier segment k

    uint NGEN;  // number of generic boundaries
    uint NNGN;  // total number of nodes for generic boundaries

    std::vector<std::vector<uint>> NBGN;  // node numbers for a segment of generic boundary

    AdcircFormat(const std::string& in_name);

    void write_to(const char* out_name) const;

    SWE::BoundaryConditions get_ibtype(std::array<uint, 2>& node_pair) const;

    bool has_edge(std::vector<uint>::const_iterator cbegin,
                  std::vector<uint>::const_iterator cend,
                  std::array<uint, 2>&              node_pair) const;
};

#endif
