#include "adcirc_format.hpp"

AdcircFormat::AdcircFormat(const std::string& fort14) {
    if (!Utilities::file_exists(fort14)) {
        throw std::logic_error("Fatal Error: ADCIRC mesh file " + fort14 + " was not found!\n");
    }

    std::ifstream ifs(fort14);

    ifs >> this->name;
    ifs.ignore(1000, '\n');

    uint n_nodes, n_elements;

    ifs >> n_elements;
    ifs >> n_nodes;
    ifs.ignore(1000, '\n');

    {  // read in node information
        uint node_name;
        std::array<double, 3> node_data;

        for (uint i = 0; i < n_nodes; ++i) {
            ifs >> node_name;
            ifs >> node_data[0] >> node_data[1] >> node_data[2];

            assert(!this->nodes.count(node_name));  // Can't define a node twice
            this->nodes.insert({node_name, node_data});
            ifs.ignore(1000, '\n');
        }
    }

    {  // read in element information
        uint element_name;
        std::array<uint, 4> element_data;

        for (uint i = 0; i < n_elements; ++i) {
            ifs >> element_name;
            ifs >> element_data[0] >> element_data[1] >> element_data[2] >> element_data[3];

            assert(!this->elements.count(element_name));  // Can't define element twice
            this->elements.insert({element_name, element_data});
            ifs.ignore(1000, '\n');
        }
    }

    {  // process open boundaries
        ifs >> this->NOPE;
        ifs.ignore(1000, '\n');
        ifs >> this->NETA;
        ifs.ignore(1000, '\n');

        this->NBDV.reserve(this->NOPE);
        uint n_nodes_bdry;
        for (uint bdry = 0; bdry < this->NOPE; ++bdry) {
            ifs >> n_nodes_bdry;
            ifs.ignore(1000, '\n');
            this->NBDV.emplace_back(n_nodes_bdry);
            for (uint n = 0; n < n_nodes_bdry; ++n) {
                ifs >> this->NBDV[bdry][n];
                ifs.ignore(1000, '\n');
            }
        }
    }

    {  // process land boundaries
        ifs >> this->NBOU;
        ifs.ignore(1000, '\n');
        ifs >> this->NVEL;
        ifs.ignore(1000, '\n');

        this->NBDV.reserve(this->NBOU);
        this->IBTYPE.resize(this->NBOU);
        uint n_nodes_bdry;
        for (uint bdry = 0; bdry < this->NBOU; ++bdry) {
            ifs >> n_nodes_bdry;
            ifs >> this->IBTYPE[bdry];
            ifs.ignore(1000, '\n');

            if (this->IBTYPE[bdry] % 10 == 0 ||  // land
                this->IBTYPE[bdry] % 10 == 1 ||  // island
                this->IBTYPE[bdry] % 10 == 2) {  // flow
                // *** //
                this->NBVV.emplace_back(n_nodes_bdry);
                for (uint n = 0; n < n_nodes_bdry; ++n) {
                    ifs >> this->NBVV[bdry][n];
                    ifs.ignore(1000, '\n');
                }
            } else if (this->IBTYPE[bdry] % 10 == 4) {  // internal barrier
                this->NBVV.emplace_back(n_nodes_bdry);
                this->IBCONN[bdry]    = std::vector<uint>(n_nodes_bdry);
                this->BARINTH[bdry]   = std::vector<double>(n_nodes_bdry);
                this->BARINCFSB[bdry] = std::vector<double>(n_nodes_bdry);
                this->BARINCFSP[bdry] = std::vector<double>(n_nodes_bdry);
                for (uint n = 0; n < n_nodes_bdry; ++n) {
                    ifs >> this->NBVV[bdry][n];
                    ifs >> this->IBCONN[bdry][n];
                    ifs >> this->BARINTH[bdry][n];
                    ifs >> this->BARINCFSB[bdry][n];
                    ifs >> this->BARINCFSP[bdry][n];
                    ifs.ignore(1000, '\n');
                }
            } else if (this->IBTYPE[bdry] == 77) {  // function
                this->NBVV.emplace_back(n_nodes_bdry);
                for (uint n = 0; n < n_nodes_bdry; ++n) {
                    ifs >> this->NBVV[bdry][n];
                    ifs.ignore(1000, '\n');
                }
            } else if (this->IBTYPE[bdry] == 88) {  // outflow
                this->NBVV.emplace_back(n_nodes_bdry);
                for (uint n = 0; n < n_nodes_bdry; ++n) {
                    ifs >> this->NBVV[bdry][n];
                    ifs.ignore(1000, '\n');
                }
            } else {
                throw std::logic_error("Fatal Error: undefined boundary type in ADCIRC mesh: " +
                                       std::to_string(this->IBTYPE[bdry]) + "!\n");
            }
        }
    }

    std::string line;
    std::getline(ifs, line);

    if (line.find_first_not_of(' ') != std::string::npos) {  // process generic boundaries if there are any, i.e. not a
                                                             // white space or end of file
        std::stringstream stream;
        stream = std::stringstream(line);

        stream >> this->NGEN;
        stream.ignore(1000, '\n');
        ifs >> this->NNGN;
        ifs.ignore(1000, '\n');

        this->NBGN.reserve(this->NGEN);
        uint n_nodes_bdry;
        for (uint bdry = 0; bdry < this->NGEN; ++bdry) {
            ifs >> n_nodes_bdry;
            ifs.ignore(1000, '\n');
            this->NBGN.emplace_back(n_nodes_bdry);
            for (uint n = 0; n < n_nodes_bdry; ++n) {
                ifs >> this->NBGN[bdry][n];
                ifs.ignore(1000, '\n');
            }
        }
    }

    ifs.close();
}

void AdcircFormat::write_to(const char* out_name) const {
    std::ofstream file;
    file.open(out_name);

    file << this->name << '\n';
    file << this->elements.size() << "    " << this->nodes.size() << '\n';

    for (const auto& node : this->nodes) {
        file << node.first;
        for (uint i = 0; i < 3; ++i) {
            file << std::setprecision(15) << std::setw(22) << node.second[i];
        }
        file << '\n';
    }

    for (const auto& element : this->elements) {
        file << element.first;
        for (uint i = 0; i < 4; ++i) {
            file << std::setw(12) << element.second[i];
        }
        file << '\n';
    }

    file << this->NOPE << " = Number of open boundaries\n";
    file << this->NETA << " = Total number of open boundary nodes\n";
    for (uint n = 0; n < this->NOPE; ++n) {
        file << this->NBDV[n].size() << " = Number of nodes for open boundary " << n + 1 << " \n";
        for (uint i : this->NBDV[n]) {
            file << i << "\n";
        }
    }

    file << this->NBOU << " = Number of land Boundaries\n";
    file << this->NVEL << " = Total number of open boundary nodes\n";
    for (uint n = 0; n < this->NBOU; ++n) {
        file << this->NBVV[n].size() << " " << this->IBTYPE[n] << " = Number of nodes for land boundary " << n + 1
             << '\n';
        for (uint i = 0; i < this->NBVV[n].size(); ++i) {
            if (this->IBTYPE[n] % 10 == 0 ||  // land
                this->IBTYPE[n] % 10 == 1 ||  // island
                this->IBTYPE[n] % 10 == 2) {  // flow
                file << this->NBVV[n][i] << '\n';
            } else if (this->IBTYPE[n] % 10 == 4) {  // internal barrier
                file << this->NBVV[n][i] << ' ' << this->IBCONN.at(n)[i] << ' ' << this->BARINTH.at(n)[i] << ' '
                     << this->BARINCFSB.at(n)[i] << ' ' << this->BARINCFSP.at(n)[i] << '\n';
            } else if (this->IBTYPE[n] == 77) {  // function
                file << this->NBVV[n][i] << '\n';
            } else if (this->IBTYPE[n] == 88) {  // outflow
                file << this->NBVV[n][i] << '\n';
            }
        }
    }

    file << this->NGEN << " = Number of generic boundaries\n";
    file << this->NNGN << " = Total number of generic boundary nodes\n";
    for (uint n = 0; n < this->NGEN; ++n) {
        file << this->NBGN[n].size() << " = Number of nodes for generic boundary " << n + 1 << " \n";
        for (uint i : this->NBGN[n]) {
            file << i << "\n";
        }
    }
}

SWE::BoundaryTypes AdcircFormat::get_ibtype(std::array<uint, 2>& node_pair) const {
    for (auto& open_bdry : this->NBDV) {
        if (has_edge(open_bdry.cbegin(), open_bdry.cend(), node_pair)) {
            return SWE::BoundaryTypes::tide;
        }
    }

    for (uint segment_id = 0; segment_id < this->NBOU; ++segment_id) {
        if (this->IBTYPE[segment_id] % 10 == 0) {
            if (has_edge(this->NBVV[segment_id].cbegin(), this->NBVV[segment_id].cend(), node_pair)) {
                return SWE::BoundaryTypes::land;
            }
        } else if (this->IBTYPE[segment_id] % 10 == 2) {
            if (has_edge(this->NBVV[segment_id].cbegin(), this->NBVV[segment_id].cend(), node_pair)) {
                return SWE::BoundaryTypes::flow;
            }
        } else if (this->IBTYPE[segment_id] % 10 == 4) {
            if (has_edge(this->NBVV[segment_id].cbegin(), this->NBVV[segment_id].cend(), node_pair) ||
                has_edge(this->IBCONN.at(segment_id).cbegin(), this->IBCONN.at(segment_id).cend(), node_pair)) {
                return SWE::BoundaryTypes::levee;
            }
        } else if (this->IBTYPE[segment_id] == 77) {
            if (has_edge(this->NBVV[segment_id].cbegin(), this->NBVV[segment_id].cend(), node_pair)) {
                return SWE::BoundaryTypes::function;
            }
        } else if (this->IBTYPE[segment_id] == 88) {
            if (has_edge(this->NBVV[segment_id].cbegin(), this->NBVV[segment_id].cend(), node_pair)) {
                return SWE::BoundaryTypes::outflow;
            }
        }
    }

    for (auto& generic_bdry : this->NBGN) {
        if (has_edge(generic_bdry.cbegin(), generic_bdry.cend(), node_pair)) {
            return SWE::BoundaryTypes::land;
        }
    }

    if (is_levee_end(node_pair)) {
        return SWE::BoundaryTypes::land;
    }

    throw std::logic_error("Fatal Error: boundary not found, unable to assign BOUNDARY_TYPE to given node_pair (" +
                           std::to_string(node_pair[0]) + ", " + std::to_string(node_pair[1]) + ")!\n");
}

std::array<uint, 2> AdcircFormat::get_internal_node_pair(std::array<uint, 2>& node_pair) const {
    for (const auto& it : this->IBCONN) {
        uint segment_id = it.first;

        auto& segment_nbvv   = this->NBVV[segment_id];
        auto& segment_ibconn = it.second;

        uint n_nodes = segment_nbvv.size();

        // Look through segment NBVV for pair
        if (segment_nbvv[0] == node_pair[0] && segment_nbvv[1] == node_pair[1]) {
            return std::array<uint, 2>{segment_ibconn[1], segment_ibconn[0]};
        }

        for (uint node = 1; node < n_nodes - 1; ++node) {
            if (segment_nbvv[node] == node_pair[0]) {
                if (segment_nbvv[node + 1] == node_pair[1]) {  // look node after
                    return std::array<uint, 2>{segment_ibconn[node + 1], segment_ibconn[node]};
                } else if (segment_nbvv[node - 1] == node_pair[1]) {  // look node before
                    return std::array<uint, 2>{segment_ibconn[node - 1], segment_ibconn[node]};
                }
            }
        }

        if (segment_nbvv[n_nodes - 1] == node_pair[0] && segment_nbvv[n_nodes - 2] == node_pair[1]) {
            return std::array<uint, 2>{segment_ibconn[n_nodes - 2], segment_ibconn[n_nodes - 1]};
        }

        // Look through segment IBCONN for pair
        if (segment_ibconn[0] == node_pair[0] && segment_ibconn[1] == node_pair[1]) {
            return std::array<uint, 2>{segment_nbvv[1], segment_nbvv[0]};
        }

        for (uint node = 1; node < n_nodes - 1; ++node) {
            if (segment_ibconn[node] == node_pair[0]) {
                if (segment_ibconn[node + 1] == node_pair[1]) {  // look node after
                    return std::array<uint, 2>{segment_nbvv[node + 1], segment_nbvv[node]};
                } else if (segment_ibconn[node - 1] == node_pair[1]) {  // look node before
                    return std::array<uint, 2>{segment_nbvv[node - 1], segment_nbvv[node]};
                }
            }
        }

        if (segment_ibconn[n_nodes - 1] == node_pair[0] && segment_ibconn[n_nodes - 2] == node_pair[1]) {
            return std::array<uint, 2>{segment_nbvv[n_nodes - 2], segment_nbvv[n_nodes - 1]};
        }
    }

    throw std::logic_error("Fatal Error: back nodes not found to given node_pair (" + std::to_string(node_pair[0]) +
                           ", " + std::to_string(node_pair[1]) + ")!\n");
}

bool AdcircFormat::has_edge(std::vector<uint>::const_iterator cbegin,
                            std::vector<uint>::const_iterator cend,
                            std::array<uint, 2>& node_pair) const {
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

bool AdcircFormat::is_levee_end(const std::array<uint, 2>& node_pair) const {
    uint min_ = std::min(node_pair[0], node_pair[1]);
    uint max_ = std::max(node_pair[0], node_pair[1]);

    for (uint segment_id = 0; segment_id < this->NBOU; ++segment_id) {
        if (this->IBTYPE[segment_id] % 10 == 4) {
            if ((min_ == std::min(*this->NBVV[segment_id].cbegin(), *this->IBCONN.at(segment_id).cbegin()) &&
                 max_ == std::max(*this->NBVV[segment_id].cbegin(), *this->IBCONN.at(segment_id).cbegin())) ||
                (min_ == std::min(*(this->NBVV[segment_id].cend() - 1), *(this->IBCONN.at(segment_id).cend() - 1)) &&
                 max_ == std::max(*(this->NBVV[segment_id].cend() - 1), *(this->IBCONN.at(segment_id).cend() - 1)))) {
                return true;
            }
        }
    }

    return false;
}