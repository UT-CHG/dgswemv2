#include "adcirc_format.hpp"

AdcircFormat::AdcircFormat(const std::string& fort14) {
    std::ifstream ifs(fort14);

    if (!ifs) {
        std::string err_msg = "Fatal Error: Mesh named " + fort14 + " not found\n";
        throw std::logic_error(err_msg);
    }

    std::getline(ifs, name);
    int n_nodes, n_elements;

    ifs >> n_elements;
    ifs >> n_nodes;
    ifs.ignore(1000, '\n');

    {  // read in node information
        int node_name;
        std::array<double, 3> node_data;

        for (int i = 0; i < n_nodes; ++i) {
            ifs >> node_name;
            ifs >> node_data[0] >> node_data[1] >> node_data[2];

            assert(!nodes.count(node_name));  // Can't define a node twice
            nodes.insert({node_name, node_data});
            ifs.ignore(1000, '\n');
        }
    }

    {  // read in element information
        int element_name;
        std::array<int, 4> element_data;

        for (int i = 0; i < n_elements; ++i) {
            ifs >> element_name;
            ifs >> element_data[0] >> element_data[1] >> element_data[2] >> element_data[3];

            assert(!elements.count(element_name));  // Can't define element twice
            elements.insert({element_name, element_data});
            ifs.ignore(1000, '\n');
        }
    }

    {  // check if element nodes are ccw, swap if necessary
        for (auto& elt : elements) {
            Shape::StraightTriangle triangle({{nodes.at(elt.second[1])[0], nodes.at(elt.second[1])[1]},
                                              {nodes.at(elt.second[2])[0], nodes.at(elt.second[2])[1]},
                                              {nodes.at(elt.second[3])[0], nodes.at(elt.second[3])[1]}});

            // note that point selection doesn't matter since the Jacobian is
            // constant
            if (!triangle.CheckJacobianPositive(Point<2>())) {
                std::swap(elt.second[1], elt.second[3]);
            }
        }
    }

    {  // process open boundaries
        ifs >> NOPE;
        ifs.ignore(1000, '\n');
        ifs >> NETA;
        ifs.ignore(1000, '\n');

        NBDV.reserve(NOPE);
        int n_nodes_bdry;
        for (int bdry = 0; bdry < NOPE; ++bdry) {
            ifs >> n_nodes_bdry;
            ifs.ignore(1000, '\n');
            NBDV.push_back(std::vector<int>(n_nodes_bdry));

            for (int n = 0; n < n_nodes_bdry; ++n) {
                ifs >> NBDV[bdry][n];
                ifs.ignore(1000, '\n');
            }
        }
    }

    {  // process land boundaries
        ifs >> NBOU;
        ifs.ignore(1000, '\n');
        ifs >> NVEL;
        ifs.ignore(1000, '\n');

        IBTYPE.resize(NBOU);
        NBVV.reserve(NBOU);
        int n_nodes_bdry;
        for (int bdry = 0; bdry < NBOU; ++bdry) {
            ifs >> n_nodes_bdry;
            ifs >> IBTYPE[bdry];
            ifs.ignore(1000, '\n');

            NBVV.push_back(std::vector<int>(n_nodes_bdry));
            for (int n = 0; n < n_nodes_bdry; ++n) {
                ifs >> NBVV[bdry][n];
                ifs.ignore(1000, '\n');
            }
        }
    }

    ifs.close();
}

void AdcircFormat::write_to(const char* out_name) const {
    std::ofstream file;
    file.open(out_name);

    file << name << '\n';
    file << elements.size() << "    " << nodes.size() << '\n';

    for (const auto& node : nodes) {
        file << node.first;
        for (int i = 0; i < 3; ++i) {
            file << std::setprecision(15) << std::setw(22) << node.second[i];
        }
        file << '\n';
    }

    for (const auto& element : elements) {
        file << element.first;
        for (int i = 0; i < 4; ++i) {
            file << std::setw(12) << element.second[i];
        }
        file << '\n';
    }

    file << NOPE << " = Number of open boundaries\n";
    file << NETA << " = Total number of open boundary nodes\n";
    for (int n = 0; n < NOPE; ++n) {
        file << NBDV[n].size() << " = Number of nodes for open boundry " << n + 1 << " \n";
        for (uint i = 0; i < NBDV[n].size(); ++i) {
            file << NBDV[n][i] << "\n";
        }
    }

    file << NBOU << " = Number of land Boundaries\n";
    file << NVEL << " = Total number of open boundary nodes\n";

    for (int n = 0; n < NBOU; ++n) {
        file << NBVV[n].size() << " " << IBTYPE[n] << " = Number of nodes for land boundary " << n + 1 << '\n';
        for (uint i = 0; i < NBVV[n].size(); ++i) {
            file << NBVV[n][i] << '\n';
        }
    }
}

SWE::BoundaryConditions AdcircFormat::get_ibtype(std::array<int, 2>& node_pair) const {
    for (auto& open_bdry : NBDV) {
        if (has_edge(open_bdry.cbegin(), open_bdry.cend(), node_pair)) {
            return SWE::BoundaryConditions::tidal;
        }
    }

    uint segment_id = 0;
    for (auto& land_bdry : NBVV) {
        if (has_edge(land_bdry.cbegin(), land_bdry.cend(), node_pair)) {
            switch (IBTYPE[segment_id]) {
                case 0:
                case 10:
                case 20:
                    return SWE::BoundaryConditions::land;
                case 2:
                case 12:
                case 22:
                    return SWE::BoundaryConditions::flow;
                default:
                    throw std::logic_error(
                        "Error Boundary type unknown, unable to assign BOUNDARY_TYPE to given "
                        "node_pair (" +
                        std::to_string(node_pair[0]) + ", " + std::to_string(node_pair[1]) + ")\n");
            }
        }

        segment_id++;
    }

    throw std::logic_error(
        "Error Boundary not found, unable to assign BOUNDARY_TYPE to given "
        "node_pair (" +
        std::to_string(node_pair[0]) + ", " + std::to_string(node_pair[1]) + ")\n");
}

bool AdcircFormat::has_edge(std::vector<int>::const_iterator cbegin,
                            std::vector<int>::const_iterator cend,
                            std::array<int, 2>& node_pair) const {
    auto it = std::find(cbegin, cend, node_pair[0]);

    if (it != cend) {
        // look at next element
        if ((it + 1) != cend && *(it + 1) == node_pair[1]) {
            return true;
        }

        // look at previous element unless we are at the first element
        if (it != cbegin) {
            if (*(it - 1) == node_pair[1]) {
                return true;
            }
        }

        // here you're not at the end of vector, but haven't found an edge
        // so we must look through the tail
        return has_edge(it + 1, cend, node_pair);
    }

    return false;
}
