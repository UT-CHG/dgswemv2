#ifndef ADCIRC_FORMAT_HPP
#define ADCIRC_FORMAT_HPP

#include "../general_definitions.hpp"
#include "../problem/SWE/swe_definitions.hpp"

class AdcircFormat {
public:
    AdcircFormat(const char* in_name);

    std::string name;
    std::unordered_map<int, std::array<double, 3>> nodes;
    std::unordered_map<int, std::array<int, 4>> elements;

    //see http://adcirc.org/home/documentation/users-manual-v51/input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14/
    int NOPE; //number of open boundaries
    int NETA; //total number of nodes for open boundary

    std::vector<std::vector<int> > NBDV; //Ordered node numbers

    int NBOU;
    int NVEL;

    std::vector<int> IBTYPE; //boundary type
    std::vector<std::vector<int> > NBVV; //node numbers on normal flow boundary segment k

public:
    SWE::BoundaryConditions get_ibtype(std::array<int, 2>& node_pair) const {
        for (auto& open_bdry : NBDV) {
            if (has_edge(open_bdry.cbegin(), open_bdry.cend(), node_pair)) {
                return SWE::BoundaryConditions::tidal;
            }
        }

        for (auto& land_bdry : NBVV) {
            if (has_edge(land_bdry.cbegin(), land_bdry.cend(), node_pair)) {
                return SWE::BoundaryConditions::land;
            }
        }

        throw std::logic_error("Error Boundary not found, unable to assign BOUNDARY_TYPE to given node_pair (" + std::to_string(node_pair[0]) + ", " + std::to_string(node_pair[1]) + ")\n");
    }

    void write_to(const char* out_name);

private:
    bool has_edge(std::vector<int>::const_iterator cbegin, std::vector<int>::const_iterator cend,
        std::array<int, 2>& node_pair) const
    {
        auto it = std::find(cbegin, cend, node_pair[0]);
        if (it != cend) {
            ++it; //look at next element
            if (it == cend) {
                return false; // if you found the last element in that boundary
            }
            else if (*it == node_pair[1]) {
                return true; //found the edge
            }
            else {
                //here you're not at the end of vector, but haven't found an edge
                //so we must look through the tail
                return has_edge(it, cend, node_pair);
            }
        }

        return false;
    }
};

#endif