#include "adcirc_format.hpp"

AdcircFormat::AdcircFormat(const std::string& fort14) {
    std::ifstream ifs(fort14);

    if (!ifs) {
        std::string err_msg = "Fatal Error: Mesh named " + fort14 + " not found\n";
        throw std::logic_error(err_msg);
    }

    std::getline(ifs, this->name);
    uint n_nodes, n_elements;

    ifs >> n_elements;
    ifs >> n_nodes;
    ifs.ignore(1000, '\n');

    {  // read in node information
        uint                  node_name;
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
        uint                element_name;
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
            this->NBDV.push_back(std::vector<uint>(n_nodes_bdry));

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

        this->IBTYPE.resize(this->NBOU);
        this->NBVV.reserve(this->NBOU);
        uint n_nodes_bdry;
        for (uint bdry = 0; bdry < this->NBOU; ++bdry) {
            ifs >> n_nodes_bdry;
            ifs >> this->IBTYPE[bdry];
            ifs.ignore(1000, '\n');

            this->NBVV.push_back(std::vector<uint>(n_nodes_bdry));
            for (uint n = 0; n < n_nodes_bdry; ++n) {
                ifs >> this->NBVV[bdry][n];
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
        file << this->NBDV[n].size() << " = Number of nodes for open boundry " << n + 1 << " \n";
        for (uint i = 0; i < this->NBDV[n].size(); ++i) {
            file << this->NBDV[n][i] << "\n";
        }
    }

    file << this->NBOU << " = Number of land Boundaries\n";
    file << this->NVEL << " = Total number of open boundary nodes\n";

    for (uint n = 0; n < this->NBOU; ++n) {
        file << this->NBVV[n].size() << " " << this->IBTYPE[n] << " = Number of nodes for land boundary " << n + 1
             << '\n';
        for (uint i = 0; i < this->NBVV[n].size(); ++i) {
            file << this->NBVV[n][i] << '\n';
        }
    }
}

SWE::BoundaryConditions AdcircFormat::get_ibtype(std::array<uint, 2>& node_pair) const {
    for (auto& open_bdry : this->NBDV) {
        if (has_edge(open_bdry.cbegin(), open_bdry.cend(), node_pair)) {
            return SWE::BoundaryConditions::tidal;
        }
    }

    uint segment_id = 0;
    for (auto& land_bdry : this->NBVV) {
        if (has_edge(land_bdry.cbegin(), land_bdry.cend(), node_pair)) {
            switch (this->IBTYPE[segment_id]) {
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

bool AdcircFormat::has_edge(std::vector<uint>::const_iterator cbegin,
                            std::vector<uint>::const_iterator cend,
                            std::array<uint, 2>&              node_pair) const {
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
